#' Wrapper function returning Summarized Experiment object giving promoter 
#' counts and activity
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
#'   and activity with gene expression stored as metadata. rowData contains
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
#' @import SummarizedExperiment
#' @importFrom data.table as.data.table .N ':='
#' @importFrom rlang .data
#' @importFrom S4Vectors metadata
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
                                                    fileLabels, drop = FALSE]))
    metadata(result) <- list(geneExpression = 
                                    geneExpression[, fileLabels, drop = FALSE])
    print('Calculating positions of promoters...')
    promoterCoordinates <- promoterCoordinates(promoterAnnotation)
    promoterIdMapping <- promoterIdMapping(promoterAnnotation)
    promoterCoordinates$geneId <- promoterIdMapping$geneId[match(
                promoterCoordinates$promoterId, promoterIdMapping$promoterId)]
    promoterCoordinates <- as.data.table(promoterCoordinates)
    promoterPosition <- NULL
    geneId <- NULL
    promoterCoordinates[, promoterPosition := ifelse(strand == '+', seq_len(.N), 
                                                rev(seq_len(.N))), by=geneId]
    ## Build row data
    rowData(result) <- data.frame(
                        absolutePromoterActivity[,c('promoterId', 'geneId')], 
                        promoterCoordinates[,c("seqnames","start", "strand",
                                    "internalPromoter", "promoterPosition")])
    transcriptByGene <- split(promoterIdMapping$transcriptName, 
                                promoterIdMapping$geneId)
    rowData(result)$txId <- transcriptByGene[match(rowData(result)$geneId, 
                                                names(transcriptByGene))]
    return(result)
}


# Helper function to summarize results across condition
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SummarizedExperiment rowData colData assays
summarizeAcrossCondition <- function(result, condition) {
        colData(result) <- DataFrame(sampleName = colnames(result), 
                                    condition=condition)
        colnames(result) <- colData(result)$sampleName
        print('Summarising expression and activity across conditions...')
        for (group in unique(condition)) {
            rowData(result)[,paste0(group, '.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$abs) 
            metadata(result)$geneExpression[,paste0(group, '.mean')] <- 
                rowMeans(metadata(result)$geneExpression[
                    ,colData(result)$condition==group, drop=FALSE])
        }
        rowData(result) <- categorizePromoters(rowData(result), condition)
    return(result)
}

# Helper function to categorize promoters
#' @importFrom data.table as.data.table .I setorder
categorizePromoters <- function(rdata, condition) {
    ## Avoid undefined global vars
    promoterId <- NULL
    promoterPosition <- NULL
    geneId <- NULL
    data <- as.data.table(rdata)
    data <- setorder(data, promoterPosition)
    data <- setorder(data, geneId)
    
    for (group in unique(condition)) {
        print(paste0('Categorizing ', group, ' promoters...'))
        mean <- paste0(group, '.mean')
        class <- paste0(group, '.class')
        
        max.rows <- data[, .I[which.max(get(mean))], by=geneId]
        data[[class]] <- ifelse(data[[mean]] < 0.25, 'Inactive', 'Minor')
        data[[class]][max.rows$V1] <- 'Major'
        data[[class]][which(data[[mean]] < 0.25)] <- 'Inactive'
        data[[class]][which(data$internalPromoter)] <- NA
    }
    data <- setorder(data, promoterId)
    return(data)
}
