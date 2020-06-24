#' Calculate the total number of junction reads overlapping with the introns of
#' each promoter for the input junction file
#'
#' @param promoterCoordinates A GRanges object containing reduced exon ranges by gene
#' @param intronRanges A Granges object containing the annotated unique intron
#'   ranges. These ranges will be used for counting the reads
#' @param file character path for the input junction bed or bam file
#' @param fileType character type of the junction bed file. Either 'tophat',
#'   'star' or 'bam'
#' @param genome character genome version
#'
#' @export
#' @return The total number of junction reads overlapping with each promoter for
#'   the input annotated intron ranges
#'
#' @examples
#' \dontrun{
#' file <- './sample1-tophat.bed'
#' junctionCounts <- calculateJunctionReadCounts(promoterCoordinates,
#'                                               intronRanges,
#'                                               file,
#'                                               fileType = 'tophat')
#' file <- './sample1.bam'                                                
#' junctionCounts <- calculateJunctionReadCounts(promoterCoordinates,
#'                                               intronRanges,
#'                                               file,
#'                                               fileType = 'bam',
#'                                               genome = 'hg19')}
#'
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#'
calculateJunctionReadCounts <- function(promoterCoordinates, intronRanges, file = '', fileType = '', genome = '') {
  
  if(fileType == 'tophat') {
    print(paste0('Processing: ', file))
    junctionTable <- readTopHatJunctions(file)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    print('File loaded into memory')
  } else if(fileType == 'star') {
    print(paste0('Processing: ', file))
    junctionTable <- readSTARJunctions(file)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    junctionTable$score <- junctionTable$um_reads  # to match the tophat style, uniquely mapped reads are used as score
    print('File loaded into memory')
  } else if (fileType == 'bam') {
    print(paste0('Processing: ', file))
    rawBam <- readGAlignments(file)
    bam <- keepStandardChromosomes(rawBam, pruning.mode = 'coarse')
    rm(rawBam)
    gc()
    junctions <- summarizeJunctions(bam, genome = genome)
    rm(bam)
    gc()
    junctionTable <- keepStandardChromosomes(junctions, pruning.mode = 'coarse')
    strand(junctionTable) <- junctionTable$intron_strand
    junctionTable <- junctionTable[,c('score')]
    print('Junctions extracted from BAM file')
  }

  print('Calculating junction counts')
  intronRanges.overlap <- findOverlaps(intronRanges, junctionTable, type = 'equal')
  intronRanges$junctionCounts <- rep(0, length(intronRanges))
  intronRanges$junctionCounts[queryHits(intronRanges.overlap)] <- junctionTable$score[subjectHits(intronRanges.overlap)]
  intronIdByPromoter <- as.vector(promoterCoordinates$intronId)
  intronId.unlist <- unlist(intronIdByPromoter)
  levels.tmp <- unique(promoterCoordinates$promoterId)
  levels.tmp <- levels.tmp[order(as.numeric(gsub('prmtr.', '', levels.tmp)))]
  promoterId.unlist <- factor(rep(promoterCoordinates$promoterId, sapply(intronIdByPromoter, length)), levels = levels.tmp)

  junctionCounts <- tapply(intronRanges$junctionCounts[intronId.unlist], promoterId.unlist, sum)
  names(junctionCounts) <- promoterCoordinates$promoterId
  return(junctionCounts)
}

#' Calculate the promoter read counts using junction read counts approach for
#' all the input junction files
#'
#' @param promoterAnnotation A PromoterAnnotation object containing
#'   intron ranges, promoter coordinates and the promoter id mapping
#' @param files A character vector. The list of junction or BAM files 
#'   for which the junction read counts will be calculated
#' @param fileLabels A character vector. The labels of junction or BAM 
#'   files for which the junction read counts will be calculated. These labels 
#'   will be used as column names for the output data.frame object
#' @param fileType A character. Type of the junction bed or bam file, either
#'   'tophat', 'star' or 'bam
#' @param genome A character. Genome version used. Must be specified if input is
#'  a BAM file. Defaults to NULL
#' @param ncores A numeric value. The number of cores to be used for
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter
#'   will be used as an argument to BiocParallel::bplapply
#'
#' @return A data.frame object. The number of junction reads per promoter (rows)
#'   for each sample (cols)
#' @export
#'
#' @examples
#' \dontrun{
#' files <- c('./sample1-tophat.bed', './sample2-tophat.bed')
#' fileLabels <- c('sample1', 'sample2')
#' promoterReadCounts <- calculatePromoterReadCounts(promoterAnnotation,
#'                                                    files,
#'                                                    fileLabels,
#'                                                    fileType = 'tophat',
#'                                                    genome = NULL
#'                                                    ncores = 1)
#'                                                    
#' files <- c('./sample1.bam', './sample2.bam')
#' fileLabels <- c('sample1', 'sample2')
#' promoterReadCounts <- calculatePromoterReadCounts(promoterAnnotation,
#'                                                    files,
#'                                                    fileLabels,
#'                                                    fileType = 'bam',
#'                                                    genome = 'hg19',
#'                                                    ncores = 1)}
#'
calculatePromoterReadCounts <- function(promoterAnnotation, files = NULL, fileLabels = NULL,
                                        fileType = NULL , genome = NULL, ncores = 1) {
  
  if (ncores > 1 & requireNamespace('BiocParallel', quietly = TRUE) == TRUE) {
    bpParameters <- BiocParallel::bpparam()
    bpParameters$workers <- ncores
    promoterReadCounts <- BiocParallel::bplapply(files, calculateJunctionReadCounts, promoterCoordinates = promoterCoordinates(promoterAnnotation),
                                                 intronRanges = intronRanges(promoterAnnotation), fileType = fileType, genome = genome, 
                                                 BPPARAM = bpParameters)
  } else {
    if (requireNamespace('BiocParallel', quietly = TRUE) == FALSE) {
      print('Warning: "BiocParallel" package is not available! Using sequential version instead...')
    }
    promoterReadCounts <- lapply(files, calculateJunctionReadCounts, promoterCoordinates = promoterCoordinates(promoterAnnotation),
                                 intronRanges = intronRanges(promoterAnnotation), fileType = fileType, genome = genome)
  }

  if (length(files) == 1) {
    promoterReadCounts <- data.frame(counts = promoterReadCounts)
  } else {
    promoterReadCounts <- as.data.frame(promoterReadCounts)
  }
  rownames(promoterReadCounts) <- promoterCoordinates(promoterAnnotation)$promoterId
  colnames(promoterReadCounts) <- fileLabels
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
