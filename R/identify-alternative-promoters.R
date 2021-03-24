#' Identifies alternative promoters.
#'
#' @param result A SummarizedExperiment object with assays giving promoter 
#'   counts and activity with gene expression stored as metadata (output from 
#'   proActiv). rowData contains promoter metadata and absolute promoter 
#'   activity summarized across conditions. Condition must be provided.
#' @param currentCondition A character vector. The current condition to be 
#'   compared. Samples corresponding to all other conditions will be compared 
#'   to this samples in this current condition.  
#' @param thresAbs A numeric value. Minimum value for promoter to be active in 
#'   absolute terms. Defaults to 0.25.
#' @param thresRel A numeric value. Minimum value for promoter to be active in 
#'   relative terms. Defaults to 0.05.
#' @param thresPval A numeric value. Adjusted p-value threshold for detecting 
#'   alternative promoters. Defaults to 0.05.
#' @param promoterFC A numeric value. Minimum fold change for a promoter in the
#'   current condition compared to all other conditions. Promoters must have at 
#'   least this magnitude of fold change for alternative usage.
#' @param geneFC A numeric value. Maximum fold change for gene expression. To
#'   identify alternative promoter usage independent of changes in gene 
#'   expression, limit the gene expression fold change.
#'
#' @export
#' @return A list of length 2. Each entry is a dataframe summarizing 
#'  up-regulated and down-regulated promoters and their corresponding genes, 
#'  if any.
#'  
#' @examples  
#' files <- list.files(system.file('extdata/vignette/junctions', 
#'                        package = 'proActiv'), 
#'                        full.names = TRUE, pattern = 'replicate5')
#' promoterAnnotation <- promoterAnnotation.gencode.v34.subset
#' result <- proActiv(files = files,
#'                        promoterAnnotation  = promoterAnnotation,
#'                        condition = rep(c('A549', 'HepG2'), each=1),
#'                        fileLabels = NULL,
#'                        ncores = 1)
#' alternativePromoters <- getAlternativePromoters(result, "A549")
#' 
#' @importFrom SummarizedExperiment rowData
getAlternativePromoters <- function(result, currentCondition,
                                    thresAbs = 0.25, thresRel = 0.05,
                                    thresPval = 0.05,
                                    promoterFC = 2.0, geneFC = 1.5) {
    
    condition <- result$condition
    if (is.null(condition)) {
        stop("The input summarized experiment must contain sample condition. 
                Run proActiv with a condition vector.")
    }
    id <- which(condition == currentCondition)
    if (length(id) == 0 ) {
        stop(paste0("Invalid input condition. Should correspond to one of: ",
                    paste(unique(condition), collapse = ' ')))
    }
    condition[id] <- currentCondition
    condition[-id] <- 'other'
    condition <- relevel(as.factor(condition), ref = 'other')
    resultAbs <- fitPromoters(result, currentCondition, 
                            type = "absolute", thres = thresAbs)
    resultRel <- fitPromoters(result, currentCondition,
                            type = "relative", thres = thresRel)
    
    altPro <- rep(0, nrow(result))
    altPro[
        resultAbs$padj < thresPval & resultRel$padj < thresPval &
        resultAbs$abs.cond > (promoterFC * resultAbs$abs.other) & 
        resultAbs$abs.cond > thresAbs &
        resultAbs$gexp.cond < (geneFC * resultAbs$gexp.other) & 
        resultAbs$gexp.cond > (resultAbs$gexp.other / geneFC)] <- 1
    altPro[
        resultAbs$padj < thresPval & resultRel$padj < thresPval &
        resultAbs$abs.other > (promoterFC * resultAbs$abs.cond) & 
        resultAbs$abs.other > thresAbs &
        resultAbs$gexp.other < (geneFC * resultAbs$gexp.cond) & 
        resultAbs$gexp.other > (resultAbs$gexp.cond / geneFC)] <- -1
    
    rdata <- data.frame(rowData(result))
    upReg <- which(altPro > 0)
    upSet <- rdata[upReg,seq_len(2)] 
    downReg <- which(altPro < 0)
    downSet <- rdata[downReg,seq_len(2)]
    
    if (length(upReg) == 0 && length(downReg) == 0) {
        message("No alternative promoters detected with current parameters. 
                Consider relaxing thresholds.")
    }
    return(list(upReg = upSet, downReg = downSet))
}

# Fits promoter activity to condition
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom S4Vectors metadata
#' @importFrom stats lm p.adjust relevel
fitPromoters <- function(result, currentCondition, type, thres) {
    if (type == "absolute") {
        raw.assay <- assays(result)$absolutePromoterActivity
    } else if (type == "relative") {
        raw.assay <- assays(result)$relativePromoterActivity
    }
    rdata <- rowData(result)
    condition <- result$condition
    nonInternalId <- which(rdata$internalPromoter == FALSE)
    pval <- rep(NaN, nrow(rdata))
    
    ## Regression - restrict to non-internal
    message(paste0("Fitting ", type, " promoter activity to condition..."))
    assay <- raw.assay[nonInternalId, ]
    num.pros <- nrow(assay)
    pval[nonInternalId] <- unlist(lapply(seq_len(num.pros), function(i) 
        tryCatch(summary(lm(unlist(assay[i,]) ~ condition))$coef[2,4], 
                warning = function(cond) 1,
                error = function(cond) NaN)))
    padj <- p.adjust(pval, method = "BH")
    
    ## Other metrics
    ## - mean absolute promoter activity each in group
    ## - gene expression
    id <- which(condition == currentCondition)
    
    mean.cond <- rowMeans(raw.assay[,id, drop = FALSE], na.rm = TRUE) 
    mean.other <- rowMeans(raw.assay[,-id, drop = FALSE], na.rm = TRUE) 
    
    gexp <- metadata(result)[[1]]
    gexp <- gexp[, result$sampleName]
    mean.gexp.cond <- rowMeans(gexp[,id, drop = FALSE], na.rm=TRUE)
    mean.gexp.cond <- mean.gexp.cond[match(rdata$geneId, rownames(gexp))]
    mean.gexp.other <- rowMeans(gexp[,-id, drop = FALSE], na.rm=TRUE)
    mean.gexp.other <- mean.gexp.other[match(rdata$geneId, rownames(gexp))]
    result <- data.frame(promoterId = rdata$promoterId,
                                    geneId = rdata$geneId,
                                    pval = pval,
                                    padj = padj,
                                    abs.cond = mean.cond,
                                    abs.other = mean.other,
                                    gexp.cond = mean.gexp.cond,
                                    gexp.other = mean.gexp.other)
    return(result)
}

