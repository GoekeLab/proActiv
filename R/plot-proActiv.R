#' Visualizes promoter activity and gene expression with boxplots
#'
#' @param result A SummarizedExperiment object return by proActiv, with assays 
#'   giving promoter counts and activity with gene expression stored as 
#'   metadata. rowData contains promoter metadata and absolute promoter 
#'   activity summarized across conditions. Condition must be provided. 
#' @param geneId A character vector. A single gene id. This identifier must
#'   correspond to the identifier in the promoter annotation. 
#' @param geneName A character vector. Common gene name to be displayed
#'   on plot. Optional. Defaults to NULL.
#' @param filterInternal A boolean. Determines if internal promoters should be 
#'   removed from the plot. Defaults to TRUE. 
#' @param col A character vector of colours to be used for plotting. 
#' 
#' @export
#' @return A list of length 3. Each entry is a plot corresponding to 
#'  absolute promoter activity, relative promoter activity and gene expression.
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
#' plots <- boxplotPromoters(result, "ENSG00000076864.19")
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata 
boxplotPromoters <- function(result, 
                                geneId, 
                                geneName = NULL, 
                                filterInternal = TRUE,
                                col = NULL) {
    
    rdata <- rowData(result)
    gexp <- as.matrix(assays(result)$gene)
    gexp <- gexp[, result$sampleName]
    rownames(gexp) <- rdata$geneId
    gexp <- gexp[!duplicated(gexp), ]
    gexp <- data.frame(gexp)
    genelst <- rdata$geneId
    
    geneIdx <- grep(geneId, genelst)
    if (length(geneIdx) == 0) {
        stop("Gene not found")
    }
    result <- result[geneIdx,]
    internalId <- TRUE
    if (filterInternal) {
        rdata <- rowData(result)
        internalId <- rdata$internalPromoter == FALSE
    }
    assay.abs <- assays(result)$abs[internalId,]
    assay.rel <- assays(result)$rel[internalId,]
    nonzero <- apply(assay.abs, 1, function(x) !all(x==0))
    assay.abs <- assay.abs[nonzero,,drop=FALSE]
    assay.rel <- assay.rel[nonzero,,drop=FALSE]
    gexp <- gexp[grep(geneId, rownames(gexp)),,drop=FALSE]
    if (nrow(assay.abs) == 0) {
        stop("Gene has no expressed non-internal promoters")
    }
    
    condition <- result$condition
    
    plot.abs <- generateBoxplot(assay.abs, condition, 
                                main = paste0("Absolute Promoter Activity ", 
                                                geneName),
                                col = col)
    plot.rel <- generateBoxplot(assay.rel, condition, 
                                main = paste0("Relative Promoter Activity ", 
                                                geneName),
                                col = col)
    plot.gexp <- generateBoxplot(gexp, condition, 
                                main = paste0("Gene Expression ", geneName), 
                                promoter = FALSE,
                                col = col)
    return(list(plot.abs, plot.rel, plot.gexp))
}

# Helper function for boxplotPromoters.
#' @importFrom data.table data.table melt
#' @importFrom ggplot2 ggplot geom_boxplot theme_light theme ggtitle aes 
#'   element_text scale_fill_manual
generateBoxplot <- function(data, 
                            condition, 
                            main, 
                            promoter = TRUE,
                            col) {
    ## Reshape to long for ggplot
    colnames(data) <- condition
    data$Feature <- rownames(data)
    data <- data.table(data)
    data <- melt(data, ncol(data))
    if (promoter) {
        data$Feature <- paste0('prmtr.', data$Feature)
    }
    names(data) <- c("Prmtr", "Condition", "Expression")
    data$Condition <- factor(data$Condition)
    
    Prmtr <- Expression <- Condition <- NULL
    ## Promoter and gene expression plots
    if (promoter) {
        outPlot <- ggplot(data = data, aes(x = Prmtr, y = Expression, fill = Condition)) + 
            geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) + 
            theme_light() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
            ggtitle(main)
    } else {
        outPlot <- ggplot(data = data, aes(x = Prmtr, y = Expression, fill = Condition)) + 
            geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) + 
            theme_light() +
            ggtitle(main)
    }
    if (!is.null(col)) {
        outPlot <- outPlot + 
            scale_fill_manual(values = col[seq_len(length(unique(condition)))])
    }
    return(outPlot)
}



#' Performs principal component analysis 
#'
#' @param result A SummarizedExperiment object return by proActiv, with assays 
#'   giving promoter counts and activity with gene expression stored as 
#'   metadata. rowData contains promoter metadata and absolute promoter 
#'   activity summarized across conditions. Condition must be provided. 
#' @param by A character vector. The assay to perform principal component 
#'  analysis by. One of promoterCounts, normalizedPromoterCounts, 
#'  absolutePromoterActivity and geneExpression (unambiguous substrings can be 
#'  supplied). Defaults to absolutePromoterActivity.  
#' @param main A character vector. Plot title (optional). Defaults to NULL.
#' @param col A vector of colours. If NULL, uses standard ggplot colours. 
#'  Defaults to NULL.
#' @param alpha A numeric value in between 0 and 1. Determines point 
#'  transparency.
#' @param cex.size A numeric value. Determines point size. 
#' 
#' @export
#' @return PCA plot.
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
#' result <- result[complete.cases(assays(result)[[1]]),]
#' plotPCA(result)
#' 
#' @importFrom ggplot2 ggplot geom_point theme_light ggtitle scale_color_manual
#'   xlab ylab xlim ylim
#' @importFrom stats prcomp
#' @importFrom SummarizedExperiment assays
plotPCA <- function(result, by = "absolutePromoterActivity",
                    main = NULL, col = NULL,
                    alpha = 0.75, cex.size = 2) {
    all <- c("promoterCounts", 
             "normalizedPromoterCounts",
             "absolutePromoterActivity",
             "geneExpression")
    assayIndex <- pmatch(by, all)
    assayIndex <- ifelse(assayIndex == 4, 5, assayIndex)
    if (is.na(assayIndex)) stop("Invalid assay.")
    message(paste0("Using assay ", names(assays(result))[assayIndex]))
    assay <- as.matrix(assays(result)[[assayIndex]])
    if (assayIndex == 5) {
        ## Remove duplicates for gene expression assay
        assay <- assay[!duplicated(rowData(result)$geneId),]
    }
    pca <- prcomp(t(assay))
    vv <- pca$sdev^2
    vv <- paste0(round(vv / sum(vv) * 100,2),"%")
    vv <- paste0("PC", seq_len(length(vv)), ": ", vv)
    pdata <- data.frame(PC1 = pca$x[,1], 
                        PC2 = pca$x[,2],
                        Condition = result$condition)
    
    lims <- max(abs(range(pca$x[,1])), abs(range(pca$x[,2]))) * c(-1,1)
    PC1 <- PC2 <- Condition <- NULL
    plot <- ggplot(pdata, aes(x = PC1, y = PC2)) + 
        geom_point(aes(color = Condition), alpha = alpha, size = cex.size) + 
        xlab(vv[1]) + ylab(vv[2]) + theme_light() + xlim(lims) + ylim(lims) 
    if (!is.null(main)) plot <- plot + ggtitle(main) 
    if (!is.null(col)) {
        nCond <- length(unique(result$condition))
        if (length(col) < nCond) {
            stop("Need as many colours as conditions.")
        }
        plot <- plot + scale_color_manual(values = col[seq_len(nCond)])
    }
    plot
}

#' Visualizes heatmap of features for samples
#'
#' @param result A SummarizedExperiment object return by proActiv, with assays 
#'   giving promoter counts and activity with gene expression stored as 
#'   metadata. rowData contains promoter metadata and absolute promoter 
#'   activity summarized across conditions. Condition must be provided.
#' @param by A character vector. The assay to visualize the heatmap for. 
#'   One of promoterCounts, normalizedPromoterCounts, 
#'   absolutePromoterActivity and geneExpression (unambiguous substrings can be 
#'   supplied). Defaults to absolutePromoterActivity.  
#' @param features Features to visualize. Either a list of promoterIds or 
#'   geneIds. The features must correspond to the assay is used, i.e., if 
#'   promoter assays are used, features must be promoterIds, while if gene
#'   expression assay is used, features must be geneIds. Defaults to NULL 
#'   (visualizes all features of the assay).
#' @param cex.legend A numeric value. Legend size.
#' @param cex.row A numeric value. Row label size.
#' @param cex.col A numeric value. Column label size.
#' @param row.margin A numeric value. Row margins.
#' @param col.margin A numeric value. Column margins.
#' @param col A vector of colours. Length should correspond to number of 
#'   experimental conditions. Defaults to NULL.
#' @param breaks A numeric vector. Breaks for heatmap plotting.
#' @param palette A character vector. One of bluered, redblue, redgreen, 
#'   greenred. Defaults to bluered. 
#' 
#' @export
#' @return Displays heatmap.
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
#' result <- result[complete.cases(assays(result)[[1]]),] 
#' plotHeatmap(result)
#' 
#' @importFrom gplots heatmap.2 redblue bluered redgreen greenred
#' @importFrom graphics legend
#' @importFrom SummarizedExperiment assays
#' @importFrom scales hue_pal 
#' @importFrom stats complete.cases quantile
plotHeatmap <- function(result, by = "absolutePromoterActivity", 
                        features = NULL, 
                        cex.legend = 0.75, cex.row = NULL, cex.col = NULL,
                        row.margin = 5, col.margin = 12,
                        col = NULL, breaks = NULL, palette = "bluered") {
    all <- c("promoterCounts", "normalizedPromoterCounts",
             "absolutePromoterActivity","geneExpression")
    assayIndex <- pmatch(by, all)
    assayIndex <- ifelse(assayIndex == 4, 5, assayIndex)
    if (is.na(assayIndex)) stop("Invalid assay.")
    message(paste0("Using assay ", names(assays(result))[assayIndex]))
    assay <- as.matrix(assays(result)[[assayIndex]])
    if (assayIndex == 5) {
        ## Remove duplicates for gene expression assay
        assay <- assay[!duplicated(rowData(result)$geneId),]
        rownames(assay) <- unique(rowData(result)$geneId)
    } else {
        rownames(assay) <- rowData(result)$promoterId
    }
    if (!is.null(features)) {
        assay <- assay[rownames(assay) %in% features, ]
        if (nrow(assay) == 0) {
            stop("Features are not found in chosen assay. Ensure promoter 
                 features are used with promoter assays and gene features used 
                 with gene expression assay.")
        }
    } else {
        message("Using all features.")
    }
    if (assayIndex < 5) rownames(assay) <- paste0("prmtr. ", rownames(assay))
    ## Map colours
    condition <- result$condition
    cond <- unique(condition)
    nCond <- length(cond)
    if (!is.null(col)) {
        if (length(col) < nCond) stop("Need as many colours as conditions.")
        cols <- col[seq_len(nCond)]
    } else {
        cols <- hue_pal()(nCond)
    }
    names(cols) <- cond
    ## Filter all rows equal to zero
    assay <- assay[rowSums(assay) > 0, ]
    assay <- t(scale(t(assay)))
    assay <- assay[complete.cases(assay),]
    ## Plot
    if (is.null(breaks)) {
        limit <- as.numeric(unlist(assay))
        lb <- as.numeric(quantile(limit, 0.05))
        ub <- as.numeric(quantile(limit, 0.95))
        breaks <- seq(lb, ub, length.out = 21)
    }
    pal <- get(palette, mode = "function")
    if (is.null(cex.row)) cex.row <- 0.2 + 1/log10(nrow(assay))
    if (is.null(cex.col)) cex.col <- 0.2 + 1/log10(ncol(assay))
    heatmap.2(assay, col = pal(length(breaks)-1), breaks = breaks, 
              trace = 'none', ColSideColors = cols[result$condition],
              margins = c(row.margin, col.margin),
              cexRow = cex.row, cexCol = cex.col)
    legend("topright", legend = cond, fill = cols, border=FALSE, bty="n", 
           cex = cex.legend)
}

