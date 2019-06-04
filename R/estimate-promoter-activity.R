#' Prepare the absolute promoter activity table including the promoter and gene
#' ids
#'
#' @param junctionReadCounts Matrix of junction read counts (rows: promoters,
#'   cols: samples)
#' @param promoterAnnotationData A PromoterAnnotation object containing the
#'   reduced exon ranges, annotated intron ranges, promoter coordinates and the
#'   promoter id mapping
#' @param log2 Logical indicating whether log2 read counts should be used
#'   (default: TRUE) or not
#' @param pseudocount Number to be used for log2 as pseudocount if log2 is TRUE
#'
#' @return data.frame of absolute promoter activity with promoter and gene ids
#' @export
#'
#' @examples
#' \dontrun{
#' absolutePromoterActivity <- getAbsolutePromoterActivity(junctionReadCounts,
#'                                                          promoterAnnotationData,
#'                                                          log2 = TRUE,
#'                                                          pseudocount = 1)
#' }
#'
#' @seealso \code{\link{preparePromoterIdMapping}} for preparing the mapping
#'   between promoters and genes, \code{\link{calculatePromoterReadCounts}} and
#'   \code{\link{normalizePromoterReadCounts}} for obtaining junction read
#'   counts
#'
getAbsolutePromoterActivity <- function(junctionReadCounts, promoterAnnotationData, log2 = TRUE, pseudocount = 1) {
  print(paste0('Calculating ', ifelse(log2 == TRUE, 'log2 ', ''), 'absolute promoter activity...'))
  conversionHelper <- unique(promoterIdMapping(promoterAnnotationData)[, c('promoterId', 'geneId')])
  conversionHelper <- conversionHelper[match(rownames(junctionReadCounts), conversionHelper$promoterId), ]
  if (log2 == TRUE & is.null(pseudocount) == FALSE) {
    junctionReadCounts <- log2(junctionReadCounts + pseudocount)
    junctionReadCounts <- as.data.frame(junctionReadCounts)
  }
  absolutePromoterActivity <- cbind(conversionHelper, junctionReadCounts)
  rownames(absolutePromoterActivity) <- conversionHelper$promoterId
  return(absolutePromoterActivity)
}

#' Prepare the gene expression table including the gene ids
#'
#' @param absolutePromoterActivity data.frame of absolute promoter activity with promoter and gene ids
#'
#' @return data.frame of gene expression with gene ids
#' @export
#'
#' @examples
#' \dontrun{
#' geneExpression <- getGeneExpression(absolutePromoterActivity)
#' }
#'
getGeneExpression <- function(absolutePromoterActivity) {
  print('Calculating gene expression...')
  conversionHelper <- absolutePromoterActivity[, (1:2)]
  if (ncol(absolutePromoterActivity) == 3) {
    geneExpression <- data.frame(counts = tapply(absolutePromoterActivity[, 3], conversionHelper$geneId, sum, na.rm = TRUE))
  } else {
    geneExpression <- as.data.frame(apply(absolutePromoterActivity[, -(1:2)], 2, tapply, conversionHelper$geneId, sum, na.rm = TRUE))
  }
  colnames(geneExpression) <- colnames(absolutePromoterActivity)[-(1:2)]
  geneExpression <- cbind(geneId = rownames(geneExpression), geneExpression)
  return(geneExpression)
}

#' Prepare the relative promoter activity table including the promoter and gene
#' ids
#'
#' @param absolutePromoterActivity data.frame of absolute promoter activity with
#'   promoter and gene ids
#' @param geneExpression data.frame of gene expression with gene ids
#'
#' @return data.frame of relative promoter activity with promoter and gene ids
#' @export
#'
#' @examples
#' \dontrun{
#' relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity,
#'                                                          geneExpression)
#' }
#'
getRelativePromoterActivity <- function(absolutePromoterActivity, geneExpression) {
  print(paste0('Calculating relative promoter activity...'))
  conversionHelper <- absolutePromoterActivity[, (1:2)]
  if (is.null(geneExpression) == TRUE) {
    geneExpression <- getGeneExpression(absolutePromoterActivity)
  }
  matchedGeneExpression <- geneExpression[match(conversionHelper$geneId, geneExpression$geneId), ]
  relativePromoterActivity <- absolutePromoterActivity[, -(1:2)] / matchedGeneExpression[, -1]
  relativePromoterActivity <- cbind(conversionHelper, relativePromoterActivity)
  return(relativePromoterActivity)
}
