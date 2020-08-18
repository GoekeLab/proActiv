#' Promoter annotation data for Gencode.v34 subsetted to 
#' chr1:10,000,000-30,000,000 including all the annotation objects
#' required for promoter activity estimation
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.v34 chr1:10,000,000-30,000,000. 
#'   The object has 3 slots: \describe{
#'   \item{intronRanges}{A GRanges object of 4,523 ranges corresponding
#'   to introns, annotated with the associated transcript.} 
#'   \item{promoterIdMapping}{The id mapping between transcript names, 
#'   promoter ids and gene ids for Gencode v34.} 
#'   \item{promoterCoordinates}{A GRanges object of 1,380 ranges
#'   showing the tss coordinate for each promoter of Gencode v34 
#'   chr1:10,000,000-30,000,000, annotated with the associated gene id, 
#'   coordinate of the 3' end of the first reduced exon, and intron id.} }
#'
"promoterAnnotation.gencode.v34.subset"
