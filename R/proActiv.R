#' Wrapper function returning Summarized Experiment object giving promoter counts and activity
#'
#' @param files A character vector. The list of input files for 
#'   which the junction read counts will be calculated
#' @param promoterAnnotation A PromoterAnnotation object containing the
#'   intron ranges, promoter coordinates and the promoter id mapping
#' @param fileLabels A character vector. The labels of input files 
#'   for which the junction read counts will be calculated. These labels will be 
#'   used as column names for each output data.frame object. If not provided,
#'   filenames will be used as labels. Defaults to NULL
#' @param genome A character. Genome version. Must be specified if input file
#'   type is a BAM file. Defaults to NULL
#' @param ncores A numeric value. The number of cores to be used for 
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
#' files <- c('./sample1-tophat.bed', './sample2-tophat.bed')
#' proActivFromJunctions <- proActiv(files = files,
#'                                   promoterAnnotation = promoterAnnotation, 
#'                                  fileLabels = NULL, 
#'                                  genome = NULL, 
#'                                  ncores = 1)
#' 
#' files <- c('./sample1.bam', './sample2.bam')
#' proActivFromBAM <- proActiv(files = files,
#'                             promoterAnnotation  = promoterAnnotation,
#'                            fileLabels = NULL,
#'                            genome = 'hg19',
#'                            ncores = 1)}
#'                            
#' @import SummarizedExperiment
#' @import S4Vectors
#'
proActiv <- function(files, promoterAnnotation, fileLabels = NULL, 
                     genome = NULL, ncores = 1) {
  
  checkFile <- file.exists(files)
  if (any(!checkFile)) {
    stop(paste0('Error: Please specify valid file paths. The following file does not exist: ', files[!checkFile]))
  }
  
  if (is.null(fileLabels)) {
    fileLabels <- make.names(tools::file_path_sans_ext(basename(files), compression = TRUE), unique = TRUE)
  }
  
  ext <- unique(tools::file_ext(files))
  if (length(ext) != 1) {
    stop("Error: More than one file type detected from given file path")
  }
  
  if (ext == 'gz' | ext == 'bz2' | ext == 'xz'){
    files.tmp <- gsub(paste0('\\.', ext), '', files)
    ext <- unique(tools::file_ext(files.tmp))
  }

  if (ext == 'bam') {
    fileType <- 'bam'
    if (is.null(genome)) {
      stop('Error: Please specify genome.')
    }
  } else if (ext == 'bed') {
    fileType <- 'tophat' 
  } else if (ext == 'junctions') {
    fileType <- 'star'
  } else {
    stop('Invalid input files: Input must either be a BAM file (.bam), 
       Tophat junctions file (.bed) or STAR junctions file (.junctions)')
  }
  
  promoterCounts <- calculatePromoterReadCounts(promoterAnnotation, files, fileLabels, 
                                                fileType, genome, ncores)
  normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
  absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts, promoterAnnotation)
  geneExpression <- getGeneExpression(absolutePromoterActivity)
  relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity, geneExpression)
  
  summarizedResults <- SummarizedExperiment(assays = list(promoterCounts = promoterCounts,
                                                          normalizedPromoterCounts = normalizedPromoterCounts,
                                                          absolutePromoterActivity = absolutePromoterActivity[, fileLabels, drop = FALSE],
                                                          relativePromoterActivity = relativePromoterActivity[, fileLabels, drop = FALSE]),
                                            rowData = absolutePromoterActivity[,c('promoterId', 'geneId')])
  metadata(summarizedResults) <- list(geneExpression = geneExpression[, fileLabels, drop = FALSE])
  
  summarizedResults 
}
