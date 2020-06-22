#' Calculate the total number of junction reads overlapping with the introns of
#' each promoter for the input junction file
#'
#' @param promoterCoordinates A GRanges object containing reduced exon ranges by gene
#' @param intronRanges A Granges object containing the annotated unique intron
#'   ranges. These ranges will be used for counting the reads
#' @param junctionFilePath character path for the input junction bed or bam file
#' @param junctionType character type of the junction bed file. Either 'tophat',
#'   'star' or 'bam'
#' @param genome character genome version
#'
#' @export
#' @return The total number of junction reads overlapping with each promoter for
#'   the input annotated intron ranges
#'
#' @examples
#' \dontrun{
#' junctionFilePath <- './sample1-tophat.bed'
#' junctionCounts <- calculateJunctionReadCounts(promoterCoordinates,
#'                                                intronRanges.annotated,
#'                                                junctionFilePath,
#'                                                junctionType = 'tophat')
#' junctionFilePath <- './sample1.bam'                                                
#' junctionCounts <- calculateJunctionReadCounts(promoterCoordinates,
#'                                                intronRanges.annotated,
#'                                                junctionFilePath,
#'                                                junctionType = 'bam',
#'                                                genome = 'hg19')}
#'
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#'
calculateJunctionReadCounts <- function(promoterCoordinates, intronRanges, junctionFilePath = '', junctionType = '', genome = '') {
  
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
  } else if (junctionType == 'bam') {
    print(paste0('Processing: ', junctionFilePath))
    rawBam <- readGAlignments(junctionFilePath)
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
#' @param promoterAnnotationData A PromoterAnnotation object containing the
#'   reduced exon ranges, annotated intron ranges, promoter coordinates and the
#'   promoter id mapping
#' @param junctionFilePaths A character vector. The list of junction or BAM files 
#'   for which the junction read counts will be calculated
#' @param junctionFileLabels A character vector. The labels of junction or BAM 
#'   files for which the junction read counts will be calculated. These labels 
#'   will be used as column names for the output data.frame object
#' @param junctionType A character. Type of the junction bed or bam file, either
#'   'tophat', 'star' or 'bam
#' @param genome A character. Genome version used. Must be specified if input is
#'  a BAM file. Defaults to NULL
#' @param numberOfCores A numeric value. The number of cores to be used for
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter
#'   will be used as an argument to BiocParallel::bplapply
#'
#' @return A data.frame object. The number of junction reads per promoter (rows)
#'   for each sample (cols)
#' @export
#'
#' @examples
#' \dontrun{
#' junctionFilePaths <- c('./sample1-tophat.bed', './sample2-tophat.bed')
#' junctionFileLabels <- c('sample1', 'sample2')
#' promoterReadCounts <- calculatePromoterReadCounts(promoterAnnotationData,
#'                                                    junctionFilePaths,
#'                                                    junctionFileLabels,
#'                                                    junctionType = 'tophat',
#'                                                    genome = NULL
#'                                                    numberOfCores = 1)
#'                                                    
#' junctionFilePaths <- c('./sample1.bam', './sample2.bam')
#' junctionFileLabels <- c('sample1', 'sample2')
#' promoterReadCounts <- calculatePromoterReadCounts(promoterAnnotationData,
#'                                                    junctionFilePaths,
#'                                                    junctionFileLabels,
#'                                                    junctionType = 'bam',
#'                                                    genome = 'hg19',
#'                                                    numberOfCores = 1)}
#'
calculatePromoterReadCounts <- function(promoterAnnotationData, junctionFilePaths = NULL, junctionFileLabels = NULL,
                                        junctionType = NULL , genome = NULL, numberOfCores = 1) {
  
  if (numberOfCores > 1 & requireNamespace('BiocParallel', quietly = TRUE) == TRUE) {
    bpParameters <- BiocParallel::bpparam()
    bpParameters$workers <- numberOfCores
    promoterReadCounts <- BiocParallel::bplapply(junctionFilePaths, calculateJunctionReadCounts, promoterCoordinates = promoterCoordinates(promoterAnnotationData),
                                                 intronRanges = annotatedIntronRanges(promoterAnnotationData), junctionType = junctionType, genome = genome, 
                                                 BPPARAM = bpParameters)
  } else {
    if (requireNamespace('BiocParallel', quietly = TRUE) == FALSE) {
      print('Warning: "BiocParallel" package is not available! Using sequential version instead...')
    }
    promoterReadCounts <- lapply(junctionFilePaths, calculateJunctionReadCounts, promoterCoordinates = promoterCoordinates(promoterAnnotationData),
                                 intronRanges = annotatedIntronRanges(promoterAnnotationData), junctionType = junctionType, genome = genome)
  }

  if (length(junctionFilePaths) == 1) {
    promoterReadCounts <- data.frame(counts = promoterReadCounts)
  } else {
    promoterReadCounts <- as.data.frame(promoterReadCounts)
  }
  rownames(promoterReadCounts) <- promoterCoordinates(promoterAnnotationData)$promoterId
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
