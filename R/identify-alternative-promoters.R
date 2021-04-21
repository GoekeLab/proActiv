#' Identifies alternative promoters.
#'
#' @param result A SummarizedExperiment object with assays giving promoter 
#'   counts, activity and gene expression (output from proActiv). rowData 
#'   contains promoter metadata and absolute promoter activity summarized 
#'   across conditions. Condition must be provided.
#' @param referenceCondition A character vector. The reference condition to be 
#'   compared. Samples corresponding to all other conditions will be compared 
#'   to this samples in this current condition.  
#' @param minAbs A numeric value. Minimum value for promoter to be active in 
#'   absolute terms. Defaults to 0.25.
#' @param minRel A numeric value. Minimum value for promoter to be active in 
#'   relative terms. Defaults to 0.05.
#' @param maxPval A numeric value. Adjusted p-value threshold for detecting 
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
getAlternativePromoters <- function(result, referenceCondition,
                                    minAbs = 0.25, minRel = 0.05,
                                    maxPval = 0.05,
                                    promoterFC = 2.0, geneFC = 1.5) {
    condition <- result$condition
    if (is.null(condition)) {
        stop("The input summarized experiment must contain sample condition. 
                Run proActiv with a condition vector.")
    }
    id <- which(condition == referenceCondition)
    if (length(id) == 0 ) {
        stop(paste0("Invalid input condition. Should correspond to one of: ",
                    paste(unique(condition), collapse = ' ')))
    }
    condition[id] <- referenceCondition
    condition[-id] <- 'other'
    condition <- relevel(as.factor(condition), ref = 'other')
    resultAbs <- fitPromoters(result, referenceCondition, 
                            type = "absolute", thres = minAbs)
    resultRel <- fitPromoters(result, referenceCondition,
                            type = "relative", thres = minRel)
    
    altPro <- rep(0, nrow(result))
    altPro[
        resultAbs$padj < maxPval & resultRel$padj < maxPval &
        resultAbs$abs.cond > (promoterFC * resultAbs$abs.other) & 
        resultAbs$abs.cond > minAbs &
        resultAbs$gexp.cond < (geneFC * resultAbs$gexp.other) & 
        resultAbs$gexp.cond > (resultAbs$gexp.other / geneFC)] <- 1
    altPro[
        resultAbs$padj < maxPval & resultRel$padj < maxPval &
        resultAbs$abs.other > (promoterFC * resultAbs$abs.cond) & 
        resultAbs$abs.other > minAbs &
        resultAbs$gexp.other < (geneFC * resultAbs$gexp.cond) & 
        resultAbs$gexp.other > (resultAbs$gexp.cond / geneFC)] <- -1
    
    rdata <- data.frame(rowData(result))
    upReg <- which(altPro > 0)
    upSet <- rdata[upReg,seq_len(2)]
    upSet <- cbind(upSet, padjAbs=resultAbs[upReg,"padj"], 
                   padjRel=resultRel[upReg, "padj"])
    downReg <- which(altPro < 0)
    downSet <- rdata[downReg,seq_len(2)]
    downSet <- cbind(downSet, padjAbs=resultAbs[downReg,"padj"], 
                   padjRel=resultRel[downReg, "padj"])
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
fitPromoters <- function(result, referenceCondition, type, thres) {
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
    id <- which(condition == referenceCondition)
    
    mean.cond <- rowMeans(raw.assay[,id, drop = FALSE], na.rm = TRUE) 
    mean.other <- rowMeans(raw.assay[,-id, drop = FALSE], na.rm = TRUE) 
    
    gexp <- assays(result)$gene
    gexp <- as.matrix(gexp[, result$sampleName])
    rownames(gexp) <- rowData(result)$geneId
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

