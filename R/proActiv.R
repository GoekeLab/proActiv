#' Wrapper function returning Summarized Experiment object giving promoter counts and activity
#'
#' @param promoterAnnotationData A PromoterAnnotation object containing the
#'   reduced exon ranges, annotated intron ranges, promoter coordinates and the
#'   promoter id mapping
#' @param filePaths A character vector. The list of input files for 
#'   which the junction read counts will be calculated
#' @param fileLabels A character vector. The labels of input files 
#'   for which the junction read counts will be calculated. These labels will be 
#'   used as column names for each output data.frame object. If not provided,
#'   filenames will be used as labels. Defaults to NULL
#' @param genome A character. Genome version. Must be specified if input file
#'   type is a BAM file. Defaults to NULL
#' @param numberOfCores A numeric value. The number of cores to be used for 
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter 
#'   will be used as an argument to BiocParallel::bplapply
#' 
#'
#' @export
#' @return A SummarizedExperiment object with assays giving promoter counts and activity 
#'   with gene expression stored as column data and promoter gene id mapping stored as 
#'   row data
#' 
#' @examples
#' \dontrun{
#' filePaths <- c('./sample1-tophat.bed', './sample2-tophat.bed')
#' proActivFromJunctions <- proActiv(promoterAnnotationData = promoterAnnotationData, 
#'                                  filePaths = filePaths, 
#'                                  fileLabels = NULL, 
#'                                  genome = NULL, 
#'                                  numberOfCores = 1)
#' 
#' filePaths <- c('./sample1.bam', './sample2.bam')
#' proActivFromBAM <- proActiv(promoterAnnotation Data = promoterAnnotationData,
#'                            filePaths = filePaths,
#'                            fileLabels = NULL,
#'                            genome = 'hg19',
#'                            numberOfCores = 1)}
#'                            
#' @import SummarizedExperiment
#' @import S4Vectors
#'
proActiv <- function(promoterAnnotationData, filePaths = NULL, fileLabels = NULL, 
                     genome = NULL, numberOfCores = 1) {
  
  if (is.null(filePaths)) {
    stop(paste0('Error: Please specify valid file paths!'))
  }
  
  checkFile <- file.exists(filePaths)
  if (any(!checkFile)) {
    stop(paste0('Error: Please specify valid file paths. The following file does not exist: ', filePaths[!checkFile]))
  }
  
  if (is.null(fileLabels)) {
    fileLabels <- tools::file_path_sans_ext(basename(filePaths))
  }
  
  ext <- unique(tools::file_ext(filePaths))
  
  if (length(ext) != 1) {
    stop("Error: More than one file type detected from given file path")
  } else if (ext == 'bam') {
    fileType <- 'bam'
  } else if (ext == 'bed') {
    fileType <- 'tophat' 
  } else if (ext == 'junctions') {
    fileType <- 'star'
  } else {
    stop('Invalid input files: Input must either be a BAM file (.bam), 
       Tophat junctions file (.bed) or STAR junctions file (.junctions)')
  }
  
  promoterCounts <- calculatePromoterReadCounts(promoterAnnotationData, filePaths, fileLabels, 
                                                fileType, genome, numberOfCores)
  normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
  absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts, promoterAnnotationData)
  geneExpression <- getGeneExpression(absolutePromoterActivity)
  relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity, geneExpression)
  
  summarizedResults <- SummarizedExperiment(assays = list(promoterCounts = promoterCounts,
                                                          normalizedPromoterCounts = normalizedPromoterCounts,
                                                          absolutePromoterActivity = subset(absolutePromoterActivity, select = -c(promoterId, geneId)),
                                                          relativePromoterActivity = subset(relativePromoterActivity, select = -c(promoterId, geneId))),
                                            rowData = absolutePromoterActivity[,c('promoterId', 'geneId')])
  metadata(summarizedResults) <- list(geneExpression = subset(geneExpression, select = -c(geneId)))
  
  summarizedResults 
}
