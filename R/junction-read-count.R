#' Calculate the total number of junction reads overlapping with the introns of
#' each promoter for the input junction file
#'
#' @param exonRanges A GRanges object containing reduced exon ranges by gene
#' @param intronRanges A Granges object containing the annotated unique intron
#'   ranges. These ranges will be used for counting the reads
#' @param junctionFilePath character path for the input junction bed file
#' @param junctionType character type of the junction bed file. Either 'tophat'
#'   or 'star'
#'
#' @export
#' @return The total number of junction reads overlapping with each promoter for
#'   the input annotated intron ranges
#'
#' @examples
#' \dontrun{
#' junctionFilePath <- './sample1-tophat.bed'
#' junctionCounts <- calculateJunctionReadCounts(exonReducedRanges,
#'                                                intronRanges.annotated,
#'                                                junctionFilePath,
#'                                                junctionType = 'tophat')}
#'
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#'
calculateJunctionReadCounts <- function(exonRanges, intronRanges, junctionFilePath = '', junctionType = 'tophat') {
  if(junctionType == 'tophat') {
    print(paste0('Processing: ', junctionFilePath))
    junctionTable <- readTopHatJunctions(junctionFilePath)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    print('File loaded into memory')
  } else if(junctionType == 'star') {
    print(paste0('Processing: ', junctionFilePath))
    junctionTable <- readSTARJunctions(junctionFilePath)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    junctionTable$score <- junctionTable$um_reads  # to match the tophat style, uniquely mapped reads are used as score
    print('File loaded into memory')
  }

  print('Calculating junction counts')
  intronRanges.overlap <- findOverlaps(intronRanges, junctionTable, type = 'equal')
  intronRanges$junctionCounts <- rep(0, length(intronRanges))
  intronRanges$junctionCounts[queryHits(intronRanges.overlap)] <- junctionTable$score[subjectHits(intronRanges.overlap)]
  intronIdByPromoter <- as.vector(exonRanges$intronId)
  intronId.unlist <- unlist(intronIdByPromoter)
  levels.tmp <- unique(exonRanges$promoterId)
  levels.tmp <- levels.tmp[order(as.numeric(gsub('prmtr.', '', levels.tmp)))]
  promoterId.unlist <- factor(rep(exonRanges$promoterId, sapply(intronIdByPromoter, length)), levels = levels.tmp)

  junctionCounts <- tapply(intronRanges$junctionCounts[intronId.unlist], promoterId.unlist, sum)
  names(junctionCounts) <- exonRanges$promoterId
  return(junctionCounts)
}

#' Calculate the promoter read counts using junction read counts approach for
#' all the input junction files
#'
#' @param exonRanges A GRanges object containing reduced exon ranges by gene
#' @param intronRanges A Granges object containing the annotated unique intron
#'   ranges. These ranges willbe used for counting the reads
#' @param junctionFilePaths A character vector. The list of junction files for
#'   which the junction read counts will be calculated
#' @param junctionFileLabels A character vector. The labels of junction files
#'   for which the junction read counts will be calculated. These labels will be
#'   used as column names for the output data.frame object
#' @param junctionType A character. Type of the junction bed file, either
#'   'tophat'(default) or 'star'
#' @param numberOfCores A numeric value. The number of cores to be used for
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter
#'   will be used argument to mclapply function hence require 'parallel' package
#'   to be installed
#'
#' @return A data.frame object. The number of junction reads per promoter (rows)
#'   for each sample (cols)
#' @export
#'
#' @examples
#' \dontrun{
#' junctionFilePaths <- c('./sample1-tophat.bed', './sample2-tophat.bed')
#' junctionFileLabels <- c('sample1', 'sample2')
#' promoterReadCounts <- calculatePromoterReadCounts(exonReducedRanges,
#'                                                    intronRanges.annotated,
#'                                                    junctionFilePaths,
#'                                                    junctionFileLabels,
#'                                                    junctionType = 'tophat',
#'                                                    numberOfCores = 1)
#' }
#'
calculatePromoterReadCounts <- function(exonRanges, intronRanges, junctionFilePaths = NULL, junctionFileLabels = NULL,
                                        junctionType = 'tophat', numberOfCores = 1) {
  if (!junctionType %in% c('tophat', 'star')) {
    stop(paste0('Error: Invalid junction type: ', junctionType, '! Possible values: "tophat" or "star"'))
  }
  if (is.null(junctionFilePaths)) {
    stop(paste0('Error: Please specify valid junction file paths!'))
  }

  checkFile <- file.exists(junctionFilePaths)
  if (any(!checkFile)) {
    stop(paste0('Error: Please specify valid junction file paths. The following file does not exist: ', junctionFilePaths[!checkFile]))
  }

  if (numberOfCores == 1) {
    promoterReadCounts <- lapply(junctionFilePaths, calculateJunctionReadCounts, exonRanges = exonRanges, intronRanges = intronRanges, junctionType = junctionType)
  } else {
    if (requireNamespace('parallel', quietly = TRUE) == TRUE) {
      promoterReadCounts <- parallel::mclapply(junctionFilePaths, calculateJunctionReadCounts, exonRanges = exonRanges, intronRanges = intronRanges, junctionType = junctionType, mc.cores = numberOfCores)
    } else {
      print('Warning: "parallel" package is not available! Using sequential version instead...')
      promoterReadCounts <- lapply(junctionFilePaths, calculateJunctionReadCounts, exonRanges = exonRanges, intronRanges = intronRanges, junctionType = junctionType)
    }
  }

  if (length(junctionFilePaths) == 1) {
    promoterReadCounts <- data.frame(counts = promoterReadCounts)
  } else {
    promoterReadCounts <- as.data.frame(promoterReadCounts)
  }
  rownames(promoterReadCounts) <- exonRanges$promoterId
  colnames(promoterReadCounts) <- junctionFileLabels
  return(promoterReadCounts)
}

#' Normalize promoter read counts using DESeq2
#'
#' @param promoterReadCounts A data.frame object. The number of junction reads
#'   per promoter (rows) for each sample (cols)
#'
#' @return A data.frame object. The normalized number of junction reads per
#'   promoter (rows) for each sample (cols) using DESeq2 counts function.
#'   Requires 'DESeq2' package to be installed
#' @export
#'
#' @examples
#' \dontrun{
#' normalizedPromoterReadCounts <- normalizePromoterReadCounts(promoterReadCounts)
#' }
#'
normalizePromoterReadCounts <- function(promoterReadCounts) {
  if (ncol(promoterReadCounts) == 1) {
    return(promoterReadCounts)
  }

  activePromoters <- which(!is.na(promoterReadCounts[, 1]))
  colData <- data.frame(sampleLabels = colnames(promoterReadCounts))
  rownames(colData) <- colnames(promoterReadCounts)

  print('Calculating normalized read counts...')
  if (requireNamespace('DESeq2', quietly = TRUE) == FALSE) {
    stop('Error: DESeq2 is not installed! For normalization DESeq2 is needed, please install it.')
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = promoterReadCounts[activePromoters, ], colData = colData, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  promoterReadCounts[activePromoters, ] <- DESeq2::counts(dds, normalized = TRUE)
  return(promoterReadCounts)
}
