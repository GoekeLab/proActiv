#' Prepare the absolute promoter activity table including the promoter and gene
#' ids
#'
#' @param junctionReadCounts Matrix of junction read counts (rows: promoters,
#'   cols: samples)
#' @param promoterAnnotation A PromoterAnnotation object containing the
#'   intron ranges, promoter coordinates and the promoter id mapping
#' @param log2 Logical indicating whether log2 read counts should be used
#'   (default: TRUE) or not
#' @param pseudocount Number to be used for log2 as pseudocount if log2 is TRUE
#'
#' @return data.frame of absolute promoter activity with promoter and gene ids
#'
getAbsolutePromoterActivity <- function(junctionReadCounts, promoterAnnotation,
                                        log2 = TRUE, pseudocount = 1) {
    message(paste0('Calculating ', ifelse(log2 == TRUE, 'log2 ', ''), 
                    'absolute promoter activity...'))
    promoterIdMapping <- promoterIdMapping(promoterAnnotation)
    conversionHelper <- unique(promoterIdMapping[, c('promoterId', 'geneId')])
    conversionHelper <- conversionHelper[match(rownames(junctionReadCounts), 
                                                conversionHelper$promoterId), ]
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
#' @param absolutePromoterActivity data.frame of absolute promoter activity 
#'   with promoter and gene ids
#'
#' @return data.frame of gene expression with gene ids#'
#'
getGeneExpression <- function(absolutePromoterActivity) {
    message('Calculating gene expression...')
    conversionHelper <- absolutePromoterActivity[, c('promoterId', 'geneId')]
    if (ncol(absolutePromoterActivity) == 3) {
        geneExpression <- data.frame(counts = tapply(
                                    absolutePromoterActivity[, 3], 
                                    conversionHelper$geneId, sum, na.rm = TRUE))
    } else {
        geneExpression <- as.data.frame(apply(
                                absolutePromoterActivity[, -c(1,2)], 2, tapply, 
                                conversionHelper$geneId, sum, na.rm = TRUE))
    }
    colnames(geneExpression) <- colnames(absolutePromoterActivity)[-c(1,2)]
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
#'
getRelativePromoterActivity <- function(absolutePromoterActivity, 
                                        geneExpression) {
    message(paste0('Calculating relative promoter activity...'))
    conversionHelper <- absolutePromoterActivity[, c('promoterId', 'geneId')]
    if (is.null(geneExpression) == TRUE) {
        geneExpression <- getGeneExpression(absolutePromoterActivity)
    }
    matchedGeneExpression <- geneExpression[match(conversionHelper$geneId, 
                                                    geneExpression$geneId), ]
    relativePromoterActivity <- absolutePromoterActivity[, -c(1,2)] / 
                                                    matchedGeneExpression[, -1]
    relativePromoterActivity <- cbind(conversionHelper, 
                                                    relativePromoterActivity)
    names(relativePromoterActivity) <- names(absolutePromoterActivity)
    return(relativePromoterActivity)
}
