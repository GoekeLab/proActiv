#' Wrapper function returning Summarized Experiment object giving promoter counts and activity
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
#' @return A SummarizedExperiment object with assays giving promoter counts and activity 
#'   with gene expression stored as column data and promoter gene id mapping stored as 
#'   row data
#' 
#' @examples
#' 
#' files <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), 
#'                     full.names = TRUE)
#' proActivFromBAM <- proActiv(files = files,
#'                             promoterAnnotation  = promoterAnnotation.gencode.v34,
#'                            fileLabels = NULL,
#'                            condition = c('A549', 'MCF7'),
#'                            genome = 'hg38',
#'                            ncores = 1)
#'                            
#' @import SummarizedExperiment
#' @import S4Vectors
#' @importFrom dplyr as_tibble '%>%' group_by mutate row_number 
#' @importFrom rlang .data
proActiv <- function(files, promoterAnnotation, fileLabels = NULL, 
                    condition = NULL, genome = NULL, ncores = 1) {
    parser <- parseFile(files, fileLabels, genome)
    fileLabels <- parser$fileLabels
    fileType <- parser$fileType
    
    promoterCounts <- calculatePromoterReadCounts(promoterAnnotation, files, fileLabels, fileType, genome, ncores)
    normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)
    absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts, promoterAnnotation)
    geneExpression <- getGeneExpression(absolutePromoterActivity)
    relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity, geneExpression)

    result <- SummarizedExperiment(assays = list(promoterCounts = promoterCounts,
                                                            normalizedPromoterCounts = normalizedPromoterCounts,
                                                            absolutePromoterActivity = absolutePromoterActivity[, fileLabels, drop = FALSE],
                                                            relativePromoterActivity = relativePromoterActivity[, fileLabels, drop = FALSE]))
    metadata(result) <- list(geneExpression = geneExpression[, fileLabels, drop = FALSE])
    
    print('Calculating positions of promoters...')
    promoterCoordinates <- promoterCoordinates(promoterAnnotation)
    promoterIdMapping <- promoterIdMapping(promoterAnnotation)
    promoterCoordinates$geneId <- promoterIdMapping$geneId[match(promoterCoordinates$promoterId, promoterIdMapping$promoterId)]
    promoterCoordinates <- as_tibble(promoterCoordinates) %>%
        group_by(.data$geneId) %>%
        mutate(promoterPosition = ifelse(strand == '+', row_number(), rev(row_number()))) 
    rowData(result) <- data.frame(absolutePromoterActivity[,c('promoterId', 'geneId')], 
                                    promoterCoordinates[,c("seqnames","start", "strand", "promoterPosition")])
    
    if (!is.null(condition)) {
        if (length(condition) != length(files)) {
            warning('Condition argument is invalid. 
                    Please ensure a one-to-one correspondence between each condition and each file.
                    Returning results not summarized by condition.')
            return(result)
        }
        colData(result) <- DataFrame(sampleName = colnames(result), 
                                    condition=condition)
        colnames(result) <- colData(result)$sampleName
        print('Summarising gene expression and promoter activity across conditions...')
        for (group in unique(condition)) {
            rowData(result)[,paste0(group, '.mean')] <- 
                rowMeans(assays(result[,colData(result)$condition==group])$abs) 
            metadata(result)$geneExpression[,paste0(group, '.mean')] <- 
                rowMeans(metadata(result)$geneExpression[,colData(result)$condition==group, drop=FALSE])
        }
        rowData(result) <- categorizePromoters(rowData(result), condition)
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
    } else if (ext == 'junctions') {
        fileType <- 'star'
    } else {
        stop('Invalid input files: Input must either be a BAM file (.bam), 
            Tophat junctions file (.bed) or STAR junctions file (.junctions)')
    }
    parsed <- list(fileLabels = fileLabels, fileType = fileType) 
    return(parsed)
}

# Helper function to categorize promoters
#' @import S4Vectors
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
    }
    data <- setorder(data, promoterId)
    return(data)
}
