#' Estimates promoter counts and activity in a single command
#'
#' @param files A character vector. The list of input files for 
#'   which the junction read counts will be calculated
#' @param promoterAnnotation A PromoterAnnotation object containing the
#'   intron ranges, promoter coordinates and the promoter id mapping
#' @param fileLabels A character vector. The labels of input files 
#'   for which the junction read counts will be calculated. These labels will be 
#'   used as column names for each output data.frame object. If not provided,
#'   filenames will be used as labels. Defaults to NULL
#' @param condition A character vector. The condition to which each sample
#'   belong to. Must correspond to the order of the files. If supplied, 
#'   results are summarized by condition. Defaults to NULL
#' @param genome A character. Genome version. Must be specified if input file
#'   type is a BAM file. Defaults to NULL
#' @param ncores A numeric value. The number of cores to be used for 
#'   counting junction reads. Defaults to 1 (no parallelization). This parameter 
#'   will be used as an argument to BiocParallel::bplapply
#' 
#'
#' @export
#' @return A SummarizedExperiment object with assays giving promoter counts 
#'   and activity with gene expression. rowData contains
#'   promoter metadata and absolute promoter activity summarized across
#'   conditions (if condition is provided)
#' 
#' @examples
#' 
#' files <- list.files(system.file('extdata/vignette/junctions', 
#'                        package = 'proActiv'), 
#'                        full.names = TRUE, pattern = 'replicate5')
#' promoterAnnotation <- promoterAnnotation.gencode.v34.subset
#' result <- proActiv(files = files,
#'                        promoterAnnotation  = promoterAnnotation,
#'                        condition = rep(c('A549', 'HepG2'), each=1),
#'                        fileLabels = NULL,
#'                        ncores = 1)
#'                            
proActiv <- function(files, promoterAnnotation, fileLabels = NULL, 
                    condition = NULL, genome = NULL, ncores = 1) {
    
    parser <- parseFile(files, fileLabels, genome)
    fileLabels <- parser$fileLabels
    fileType <- parser$fileType
    
    result <- buildSummarizedExperiment(promoterAnnotation, files, fileLabels,
                                        fileType, genome, ncores)
    
    if (!is.null(condition)) {
        if (length(condition) != length(files)) {
            warning('Condition argument is invalid. 
                Please ensure a 1-1 map between each condition and each file.
                Returning results not summarized across conditions.')
            return(result)
        } else {
            result <- summarizeAcrossCondition(result, condition) 
        }
    }
    return(result) 
}

#' Integrate multiple proActiv runs 
#'
#' @param res1 A summarizedExperiment object returned by proActiv
#' @param res2 A summarizedExperiment object returned by proActiv
#' @param ... Additional summarizedExperiment objects returned by proActiv
#' @param promoterAnnotation Promoter annotation object used to create 
#'   proActiv runs
#' @param renormalize Whether to renormalize counts after merging. Defaults to 
#'   TRUE
#'
#' @export
#' @return A SummarizedExperiment object with assays giving promoter counts 
#'   and activity with gene expression. rowData contains
#'   promoter metadata and absolute promoter activity summarized across
#'   conditions (if condition is provided)
#' 
#' @examples
#' f1 <- list.files(system.file('extdata/vignette/junctions', 
#'                              package = 'proActiv'), 
#'                  full.names = TRUE, pattern = 'A549')
#' f2 <- list.files(system.file('extdata/vignette/junctions', 
#'                              package = 'proActiv'), 
#'                  full.names = TRUE, pattern = 'HepG2')
#' promoterAnnotation <- promoterAnnotation.gencode.v34.subset
#' res1 <- proActiv(files = f1, promoterAnnotation  = promoterAnnotation,
#'                  condition = rep('A549',3))
#' res2 <- proActiv(files = f2, promoterAnnotation = promoterAnnotation,
#'                  condition = rep('HepG2',3))
#' res <- integrateProactiv(res1, res2, promoterAnnotation = promoterAnnotation)
#'
#' @importFrom SummarizedExperiment assays `assays<-` cbind rowData `rowData<-`
integrateProactiv <- function(res1, res2, ..., 
                              promoterAnnotation,
                              renormalize = TRUE) {
    combined <- cbind(res1, res2, ...)
    if (renormalize) {
        promoterCounts <- assays(combined)$promoterCounts
        normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
        absolutePromoterActivity <- getAbsolutePromoterActivity(
            normalizedPromoterCounts, promoterAnnotation)
        
        geneExpression <- getGeneExpression(absolutePromoterActivity)
        relativePromoterActivity <- getRelativePromoterActivity(
            absolutePromoterActivity, geneExpression)
        
        rownames(geneExpression) <- rownames(promoterCounts)
        fileLabels <- combined$sampleName
        assays(combined) <- list(promoterCounts = promoterCounts, 
            normalizedPromoterCounts = normalizedPromoterCounts, 
            absolutePromoterActivity = absolutePromoterActivity[, fileLabels, 
                                                                drop=FALSE], 
            relativePromoterActivity = relativePromoterActivity[, fileLabels, 
                                                                drop=FALSE],
            geneExpression = geneExpression[, fileLabels, drop=FALSE])    
        
        ## Settle rowData
        rowData(combined) <- rowData(combined)[-grep("mean|class", 
                                                colnames(rowData(combined)))]
        combined <- summarizeAcrossCondition(combined, combined$condition)
    }
    return(combined)
} 



# Helper function to impute file labels and infer file type
parseFile <- function(files, fileLabels, genome) {
    checkFile <- file.exists(files)
    if (any(!checkFile)) {
        stop(paste0('Error: Please specify valid file paths. 
                    The following file does not exist: ', files[!checkFile]))
    }
    
    if (is.null(fileLabels)) {
        fileLabels <- make.names(tools::file_path_sans_ext(basename(files), 
                                                            compression = TRUE),
                                unique = TRUE)
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
    } else if (ext == 'junctions' | ext == 'tab') {
        fileType <- 'star'
    } else {
        stop('Invalid input files: Input must either be a BAM file (.bam), 
            Tophat junctions file (.bed) or 
            STAR junctions file (.junctions / .tab)')
    }
    parsed <- list(fileLabels = fileLabels, fileType = fileType) 
    return(parsed)
}

# Call functions to get activity and build summarized experiment
#' @importFrom SummarizedExperiment SummarizedExperiment 'rowData<-'
#' @importFrom data.table as.data.table .N ':='
#' @importFrom rlang .data
#' @importFrom S4Vectors 'metadata<-'
buildSummarizedExperiment <- function(promoterAnnotation, 
                                        files, fileLabels, fileType, 
                                        genome, ncores) {
    promoterCounts <- calculatePromoterReadCounts(promoterAnnotation, 
                                                    files, fileLabels, fileType, 
                                                    genome, ncores)
    normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
    absolutePromoterActivity <- getAbsolutePromoterActivity(
                                                    normalizedPromoterCounts, 
                                                    promoterAnnotation)
    geneExpression <- getGeneExpression(absolutePromoterActivity)
    relativePromoterActivity <- getRelativePromoterActivity(
                                                    absolutePromoterActivity, 
                                                    geneExpression)
    result <- SummarizedExperiment(assays = list(
            promoterCounts = promoterCounts,
            normalizedPromoterCounts = normalizedPromoterCounts,
            absolutePromoterActivity = absolutePromoterActivity[, 
                                                    fileLabels, drop = FALSE],
            relativePromoterActivity = relativePromoterActivity[, 
                                                    fileLabels, drop = FALSE],
            geneExpression = geneExpression[, fileLabels, drop = FALSE]))
    
    message('Calculating positions of promoters...')
    promoterCoordinates <- promoterCoordinates(promoterAnnotation)
    promoterIdMapping <- promoterIdMapping(promoterAnnotation)
    promoterCoordinates$geneId <- promoterIdMapping$geneId[match(
                promoterCoordinates$promoterId, promoterIdMapping$promoterId)]
    promoterCoordinates <- as.data.table(promoterCoordinates)
    promoterPosition <- geneId <- strand <- NULL
    promoterCoordinates[, promoterPosition := ifelse(strand == '+', seq_len(.N), 
                                                rev(seq_len(.N))), by=geneId]
    ## Build row data
    rowData(result) <- data.frame(
                        absolutePromoterActivity[,c('promoterId', 'geneId')], 
                        promoterCoordinates[,c("seqnames","start", "strand",
                                    "internalPromoter", "promoterPosition")])
    transcriptByPromoter <- split(promoterIdMapping$transcriptName, 
                                promoterIdMapping$promoterId)
    rowData(result)$txId <- transcriptByPromoter[match(rowData(result)$promoterId, 
                                                names(transcriptByPromoter))]
    return(result)
}


# Helper function to summarize results across condition
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SummarizedExperiment rowData colData 'colData<-' assays
summarizeAcrossCondition <- function(result, condition) {
        if (any(make.names(condition) != condition)) {
            warning("Condition is modified to be syntactically valid")
            condition <- make.names(condition)
        }
        colData(result) <- DataFrame(sampleName = colnames(result), 
                                    condition=condition)
        colnames(result) <- colData(result)$sampleName
        message('Summarising expression and activity across conditions...')
        for (group in unique(condition)) {
            rowData(result)[,paste0(group, '.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$abs) 
            rowData(result)[,paste0(group, '.gene.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$gene)
        }
        rowData(result) <- categorizePromoters(rowData(result), condition)
    return(result)
}

# Helper function to categorize promoters
#' @importFrom rlang .data
#' @importFrom dplyr as_tibble '%>%' slice_max
categorizePromoters <- function(rdata, condition) {
    rdata <- as_tibble(rdata)
    for (group in unique(condition)) {
        message(paste0('Categorizing ', group, ' promoters...'))
        mean <- paste0(group, '.mean')
        class <- paste0(group, '.class')
        
        max_rows <- rdata %>%
            group_by(.data$geneId) %>%
            slice_max(!!as.name(mean), with_ties = FALSE)
        rdata[[class]] <- ifelse(rdata[[mean]] < 0.25, 'Inactive', 'Minor')
        rdata[[class]][match(max_rows$promoterId, rdata$promoterId)] <- "Major"
        rdata[[class]][which(rdata[[mean]] < 0.25)] <- "Inactive"
        rdata[[class]][which(rdata$internalPromoter)] <- NA
        rdata[[class]][which(is.na(rdata$internalPromoter))] <- NA
    }
    return(rdata)
}
