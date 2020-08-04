#' Promoter annotation data for Gencode.v19 including all the annotation objects
#' required for promoter activity estimation
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' v19
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.v19. The object has 3 slots: \describe{
#'   \item{intronRanges}{A GRanges object of 344,651 ranges corresponding
#'   to introns, annotated with the associated transcript.} 
#'   \item{promoterIdMapping}{The id mapping between transcript names, 
#'   promoter ids and gene ids for Gencode v19.} 
#'   \item{promoterCoordinates}{A GRanges object of 113,076 ranges
#'   showing the tss coordinate for each promoter of Gencode v19,
#'   annotated with the associated gene id, coordinate of the 3' end of the first
#'   reduced exon, and intron id.} }
#'
"promoterAnnotation.gencode.v19"

#' Promoter annotation data for Gencode.v34 including all the annotation objects
#' required for promoter activity estimation
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' v34
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.v34. The object has 3 slots: \describe{
#'   \item{intronRanges}{A GRanges object of 383,654 ranges corresponding
#'   to introns, annotated with the associated transcript.} 
#'   \item{promoterIdMapping}{The id mapping between transcript names, 
#'   promoter ids and gene ids for Gencode v34.} 
#'   \item{promoterCoordinates}{A GRanges object of 122,635 ranges
#'   showing the tss coordinate for each promoter of Gencode v34,
#'   annotated with the associated gene id, coordinate of the 3' end of the first
#'   reduced exon, and intron id.} }
#'
"promoterAnnotation.gencode.v34"

#' Promoter annotation data for Gencode.vM1 including all the annotation objects
#' required for promoter activity estimation
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' vM1
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.vM1. The object has 3 slots: \describe{
#'   \item{intronRanges}{A GRanges object of 243,332 ranges corresponding
#'   to introns, annotated with the associated transcript.} 
#'   \item{promoterIdMapping}{The id mapping between transcript names, 
#'   promoter ids and gene ids for Gencode vM1.} 
#'   \item{promoterCoordinates}{A GRanges object of 60,768 ranges
#'   showing the tss coordinate for each promoter of Gencode vM1,
#'   annotated with the associated gene id, coordinate of the 3' end of the first
#'   reduced exon, and intron id.} }
#'
"promoterAnnotation.gencode.vM1"

#' Promoter annotation data for Gencode.vM25 including all the annotation objects
#' required for promoter activity estimation
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' vM25
#'
#' @format A PromoterAnnotation (S4 Class) object containing all the promoter
#'   annotation objects for Gencode.vM25. The object has 3 slots: \describe{
#'   \item{intronRanges}{A GRanges object of 285,067 ranges corresponding
#'   to introns, annotated with the associated transcript.} 
#'   \item{promoterIdMapping}{The id mapping between transcript names, 
#'   promoter ids and gene ids for Gencode vM25.} 
#'   \item{promoterCoordinates}{A GRanges object of 91,902 ranges
#'   showing the tss coordinate for each promoter of Gencode vM25,
#'   annotated with the associated gene id, coordinate of the 3' end of the first
#'   reduced exon, and intron id.} }
#'
"promoterAnnotation.gencode.vM25"
