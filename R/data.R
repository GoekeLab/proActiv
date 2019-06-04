#' Promoter annotation data for Gencode.v19 including all the annotation objects
#' required for promoter activity estimation
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' v19
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.v19. The object has 4 slots: \describe{
#'   \item{reducedExonRanges}{A GRanges object of 113,076 reduced exon ranges
#'   with promoter metadata for Gencode v19} \item{annotatedIntronRanges}{A
#'   GRanges object of 344,651 ranges with metadata. The intron ranges annotated
#'   with the promoter information} \item{promoterIdMapping}{The id mapping
#'   between transcript ids, names, TSS ids, promoter ids and gene ids for
#'   Gencode v19} \item{promoterCoordinates}{A GRanges object of 113,076 ranges
#'   showing the tss coordinate for each promoter of Gencode v19} }
#'
"promoterAnnotationData.gencode.v19"
