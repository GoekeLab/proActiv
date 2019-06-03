#' Reduced exon ranges for the first exons of each gene to determine the
#' transcript to promoter mapping for Gencode v19
#'
#' A GRanges object containing the reduced first exons ranges by gene with
#' promoter metadata
#'
#' @format A GRanges object of 113,076 reduced exon ranges with promoter
#'   metadata for Gencode v19
#'
"reducedExonRanges.gencode.v19"

#' The intron ranges annotated with the promoter information.
#'
#' A GRanges object. The intron ranges annotated with the promoter information.
#'
#' @format A GRanges object of 344,651 ranges with metadata. The intron ranges
#'   annotated with the promoter information.
#'
"annotatedIntronRanges.gencode.v19"

#' The id mapping between transcript ids, names, TSS ids, promoter ids and gene
#' ids
#'
#' A data.frame object. The id mapping between transcript ids, names, TSS ids,
#' promoter ids and gene ids
#'
#' @format A data frame with 196520 rows and 5 variables: \describe{
#'   \item{transcriptId}{the transcript ids of all transcripts for Gencode v19}
#'   \item{transcriptName}{the transcript names of all transcripts for Gencode
#'   v19} \item{tssId}{the custom tss ids of each transcript for Gencode v19}
#'   \item{promoterId}{the custom promoter ids of each transcript for Gencode
#'   v19} \item{geneId}{the gene ids of all transcripts for Gencode v19} }
#'
"promoterIdMapping.gencode.v19"

#' Transcription start site coordinates for each promoter identified using
#' Gencode v19
#'
#' A GRanges object containing the tss coordinate for each promoter for Gencode
#' v19
#'
#' @format A GRanges object of 113,076 ranges showing the tss coordinate for
#'   each promoter of Gencode v19
#'
"promoterCoordinates.gencode.v19"
