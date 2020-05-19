#' Wrapper function returning Summarized Experiment object giving promoter counts and activity
#'
#' @param promoterAnnotationData A PromoterAnnotation object containing the
#'   reduced exon ranges, annotated intron ranges, promoter coordinates and the
#'   promoter id mapping
#' @param filePaths A character vector. The list of junction files for 
#'   which the junction read counts will be calculated
#' @param fileLabels A character vector. The labels of junction files 
#'   for which the junction read counts will be calculated. These labels will be 
#'   used as column names for each output data.frame object
#' @param fileType A character. Type of input data - either a bam file or a junction
#'   file. Either 'bam', 'tophat' (default) or 'star'
#' @param genome A character. Genome version
#' @param numberOfCores A numeric value. The number of cores to be used for 
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter 
#'   will be used as an argument to parallel::mclapply or BiocParallel::bplapply
#' 
#'
#' @export
#' @return A SummarizedExperiment object with assays giving promoter counts and activity 
#'   with gene expression stored as column data and promoter gene id mapping stored as 
#'   row data
#' 
#' @examples
#' \dontrun{
#' proActivFromJunctions <- proActiv(promoterAnnotationData = promoterAnnotationData, 
#'                                  filePaths = filePaths, 
#'                                  fileLabels = NULL, 
#'                                  fileType = 'tophat', 
#'                                  genome = NULL, 
#'                                  numberOfCores = 1)
#'
#' proActivFromBAM <- proActiv(promoterAnnotation Data = promoterAnnotationData,
#'                            filePaths = filePaths,
#'                            fileLabels = NULL,
#'                            fileType = 'bam',
#'                            genome = 'hg19',
#'                            numberOfCores = 1)}
#' @import SummarizedExperiment
#' @import S4Vectors
#'
proActiv <- function(promoterAnnotationData, filePaths = NULL, fileLabels = NULL, fileType = 'tophat', 
                     genome = NULL, numberOfCores = 1) {
  promoterCounts <- calculatePromoterReadCounts(promoterAnnotationData, filePaths, fileLabels, fileType,
                                                genome, numberOfCores)
  normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
  absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts, promoterAnnotationData)
  geneExpression <- getGeneExpression(absolutePromoterActivity)
  relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity, geneExpression)
  
  summarizedResults <- SummarizedExperiment(assays = list(promoterCounts = promoterCounts,
                                            normalizedPromoterCounts = normalizedPromoterCounts,
                                            absolutePromoterActivity = subset(absolutePromoterActivity, select = -c(promoterId, geneId)),
                                            relativePromoterActivity = subset(absolutePromoterActivity, select = -c(promoterId, geneId))))
  rowData(summarizedResults) <- absolutePromoterActivity[,c('promoterId', 'geneId')]
  colData(summarizedResults) <- DataFrame(t(subset(geneExpression, select = -c(geneId))))
  summarizedResults
}
