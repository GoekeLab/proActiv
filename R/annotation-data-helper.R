#' @import IRanges
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomicAlignments
#'
NULL

# Get the transcript ranges with metadata and transcript lengths
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-'
getTranscriptRanges <- function(txdb, species = 'Homo_sapiens') {
    # Retrieve transcript ranges with detailed metadata
    transcriptRanges <- transcripts(txdb, columns = c('TXID', 'TXSTART', 'TXEND', 'GENEID', 'TXNAME'))
    names(transcriptRanges) <- transcriptRanges$TXNAME
    transcriptRanges <- keepStandardChromosomes(x = transcriptRanges, species = species, pruning.mode = "tidy")
    GenomeInfoDb::seqlevelsStyle(transcriptRanges) <- 'UCSC'
    return(transcriptRanges)
}

# Get transcript start site coordinates
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-'
getTssRanges <- function(transcriptRanges) {
    # Annotate transcription start sites (TSSs) with ids
    tssCoordinates <- promoters(transcriptRanges, upstream = 0, downstream = 1)
    seqlevelsStyle(tssCoordinates) <- 'UCSC'
    tssCoordinates.unique <- unique(tssCoordinates)
    tssCoordinates.unique$tssId <- paste0('tss.', seq_len(length(tssCoordinates.unique)))
    tssCoordinates.overlap <- findOverlaps(tssCoordinates, tssCoordinates.unique, type = 'equal')
    tssCoordinates$tssId <- tssCoordinates.unique$tssId[subjectHits(tssCoordinates.overlap)]
    return(tssCoordinates)
}

# Get exon ranges for each transcript
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-'
getExonRangesByTx <- function(txdb, species = 'Homo_sapiens') {
    exonRangesByTx <- exonsBy(txdb, by = 'tx', use.names = TRUE)
    exonRangesByTx <- keepStandardChromosomes(x = exonRangesByTx, species = species, pruning.mode = "tidy")
    seqlevelsStyle(exonRangesByTx) <- 'UCSC'
    return(exonRangesByTx)
}

# Get the first exon of each transcript
getFirstExonRanges <- function(exonRangesByTx.unlist) {
    exonRanges.firstExon <- exonRangesByTx.unlist[exonRangesByTx.unlist$exon_rank == 1]
    exonRanges.firstExon$tx_name <- names(exonRanges.firstExon)
    exonRanges.firstExon$customId <- seq_len(length(exonRanges.firstExon))
    return(exonRanges.firstExon)
}

# Reduce all first exons for each gene to identify transcripts belonging to each promoter
#' @importFrom rlang .data
#' @importFrom dplyr as_tibble mutate group_by '%>%' group_by summarise ungroup select lead arrange
getReducedExonRanges <- function(exonRanges.firstExon, exonRanges.firstExon.geneId) {
    exonRanges.firstExon$geneId <- exonRanges.firstExon.geneId
    exonReducedRanges <- as_tibble(exonRanges.firstExon) %>% 
        dplyr::arrange(.data$geneId, .data$start) %>% 
        dplyr::group_by(.data$geneId) %>%
        dplyr::mutate(exonClass = c(0, cumsum(lead(.data$start) > cummax(.data$end))[-n()])) %>%
        dplyr::group_by(.data$geneId, .data$exonClass, .data$seqnames, .data$strand) %>%
        dplyr::summarise(start = min(.data$start), 
            end = max(.data$end), 
            customId = list(.data$customId),
            .groups = 'drop') %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$seqnames, .data$start, .data$end, .data$strand, .data$customId)
    exonReducedRanges$promoterId <- paste0('prmtr.', seq_len(nrow(exonReducedRanges)))
    exonReducedRanges <- makeGRangesFromDataFrame(exonReducedRanges, keep.extra.columns = TRUE)
    names(mcols(exonReducedRanges)) <- c('revmap', 'promoterId')
    return(exonReducedRanges)
}

# Get intron ranges for each promoter
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-'
getIntronRangesByTx <- function(txdb, transcriptRanges, species = 'Homo_sapiens') {
    # Retrieve intron ranges for each promoter
    intronRangesByTx <- intronsByTranscript(txdb, use.names = TRUE)
    intronRangesByTx <- keepStandardChromosomes(x = intronRangesByTx, species = species, pruning.mode = "tidy")
    seqlevelsStyle(intronRangesByTx) <- 'UCSC'
    intronRangesByTx <- intronRangesByTx[names(transcriptRanges)]
    return(intronRangesByTx)
}

# Get unique intron ranges with custom ids
getUniqueIntronRanges <- function(intronRangesByTx.unlist) {
    # Annotate unique intron ranges  with unique ids
    intronRanges.unique <- unique(intronRangesByTx.unlist)
    intronRanges.unique <- sort(intronRanges.unique)
    intronRanges.unique$INTRONID <- seq_len(length(intronRanges.unique))
    names(intronRanges.unique) <- NULL
    return(intronRanges.unique)
}

# Get the rank of each intron within the transcript (similar to exon ranges)
#' @importFrom rlang .data
#' @importFrom dplyr as_tibble mutate group_by '%>%' row_number
getIntronRanks <- function(intronRangesByTx) {
    intronRankByTx <- as_tibble(intronRangesByTx) %>% 
        dplyr::group_by(.data$group_name) %>% 
        dplyr::mutate(rank = ifelse(strand == '+', row_number(), rev(row_number())))
    intronRank <- intronRankByTx$rank
    return(intronRank)
}

# Annotate all the intron ranges with the metadata
annotateAllIntronRanges <- function(intronRankByTx, intronRangesByTx.unlist, transcriptRanges, promoterIdMapping) {
    # Annotate all intron ranges with corresponding ids
    intronRangesByTx.unlist$INTRONRANK <- intronRankByTx
    intronRangesByTx.unlist$TXID <- transcriptRanges$TXID[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$TXSTART <- transcriptRanges$TXSTART[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$TXEND <- transcriptRanges$TXEND[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$GENEID <- as.character(transcriptRanges$GENEID)[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$TXNAME <- transcriptRanges$TXNAME[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$TxWidth <- transcriptRanges$TxWidth[match(names(intronRangesByTx.unlist), transcriptRanges$TXNAME)]
    intronRangesByTx.unlist$promoterId <- promoterIdMapping$promoterId[match(intronRangesByTx.unlist$TXID, promoterIdMapping$transcriptId)]
    names(intronRangesByTx.unlist) <- NULL
    return(intronRangesByTx.unlist)
}

# Annotate unique intorn ranges with metadata
annotateUniqueIntronRanges <- function(intronRanges.unique, intronRangesByTx.unlist) {
    # Annotate unique intron ranges with corresponding ids
    intronRanges.unique$INTRONRANK <- as(split(intronRangesByTx.unlist$INTRONRANK, intronRangesByTx.unlist$INTRONID), 'IntegerList')
    intronRanges.unique$TXID <- as(split(intronRangesByTx.unlist$TXID, intronRangesByTx.unlist$INTRONID), 'IntegerList')
    intronRanges.unique$TXSTART <- as(split(intronRangesByTx.unlist$TXSTART, intronRangesByTx.unlist$INTRONID), 'IntegerList') 
    intronRanges.unique$TXEND <- as(split(intronRangesByTx.unlist$TXEND, intronRangesByTx.unlist$INTRONID), 'IntegerList')
    intronRanges.unique$GENEID <- as(split(intronRangesByTx.unlist$GENEID, intronRangesByTx.unlist$INTRONID), 'CharacterList')
    intronRanges.unique$TXNAME <- as(split(intronRangesByTx.unlist$TXNAME, intronRangesByTx.unlist$INTRONID), 'CharacterList')
    intronRanges.unique$TxWidth <- as(split(intronRangesByTx.unlist$TxWidth, intronRangesByTx.unlist$INTRONID), 'IntegerList')
    intronRanges.unique$promoterId <- as(split(intronRangesByTx.unlist$promoterId, intronRangesByTx.unlist$INTRONID), 'CharacterList')
    return(intronRanges.unique)
}

# Prepare metadata for each unique intron range considering all of its uses across transcripts
#' @importFrom rlang .data
#' @importFrom dplyr as_tibble mutate group_by '%>%' filter distinct inner_join
getIntronTable <- function(intronRanges.unique, intronRangesByTx.unlist) {
    intronTable <- as_tibble(as.data.frame(intronRanges.unique)) 
    # Prepare metadata for each intron
    intronRangesByTx.metadata <- as_tibble(as.data.frame(mcols(intronRangesByTx.unlist)))
    intronRangesByTx.metadata <- mutate(group_by(intronRangesByTx.metadata, .data$TXNAME), IntronEndRank = max(.data$INTRONRANK) - .data$INTRONRANK + 1)
    intronRangesByTx.metadata <- mutate(group_by(intronRangesByTx.metadata, .data$INTRONID), MinIntronRank = min(.data$INTRONRANK))
    intronRangesByTx.metadata <- mutate(group_by(intronRangesByTx.metadata, .data$INTRONID), MaxIntronRank = max(.data$INTRONRANK))

    intronRangesByTx.metadata.unique <- group_by(intronRangesByTx.metadata, .data$INTRONID) %>%
        filter(.data$INTRONRANK == .data$MinIntronRank) %>%
        mutate(TxWidthMax = max(.data$TxWidth)) %>%
        dplyr::select(.data$INTRONID, .data$MinIntronRank, .data$MaxIntronRank, .data$TxWidthMax) %>%
        dplyr::rename(MinIntronRankTxWidthMax = .data$TxWidthMax) %>%
        distinct()

    intronTable <- inner_join(intronTable, intronRangesByTx.metadata.unique, by = 'INTRONID')
    return(intronTable)
}

# Get the ids of introns contributing to the transcription for each promoter
getIntronIdPerPromoter <- function(intronRangesByTx.unlist, promoterIdMapping) {
    intronRanges.firstIntron <- intronRangesByTx.unlist[intronRangesByTx.unlist$INTRONRANK == 1]
    intronIdByPromoter.firstIntron <- split(intronRanges.firstIntron$INTRONID, 
                                            factor(intronRanges.firstIntron$promoterId, levels = paste0('prmtr.', seq_len(length(unique(promoterIdMapping$promoterId))))))
    intronIdByPromoter.firstIntron <- lapply(intronIdByPromoter.firstIntron, unique)
    return(intronIdByPromoter.firstIntron)
}

# Get the metadata for each promoter regarding the use of introns
#' @importFrom rlang .data
#' @importFrom dplyr group_by mutate filter distinct n
getPromoterMetadata <- function(intronTable.firstIntron) {
    ## Extract per promoter information from intronTable
    promoterMetadata <- group_by(intronTable.firstIntron, .data$promoterId) %>%
        mutate(MinMergedIntronRank = min(.data$MinIntronRank),
            MaxMergedIntronRank = max(.data$MaxIntronRank),
            MaxMergedIntronTxWidth = max(.data$MinIntronRankTxWidthMax)) %>%
        filter(.data$MinIntronRank == min(.data$MinIntronRank)) %>%
        dplyr::select(.data$promoterId, .data$MinMergedIntronRank, .data$MaxMergedIntronRank, .data$MaxMergedIntronTxWidth) %>%
        distinct()
    return(promoterMetadata)
}

# Annotate the reduced exon ranges for each promoter with the intron metadata
annotateReducedExonRanges <- function(exonReducedRanges, promoterMetadata, intronIdByPromoter.firstIntron) {
    intronReducedTable <- cbind(as.data.frame(exonReducedRanges), as.data.frame(promoterMetadata[match(exonReducedRanges$promoterId, promoterMetadata$promoterId), -1]))
    exonReducedRanges <- GRanges(intronReducedTable$seqnames,
                                ranges = IRanges(start = intronReducedTable$start, intronReducedTable$end),
                                strand = intronReducedTable$strand)
    exonReducedRanges$intronId <- as(intronIdByPromoter.firstIntron, 'IntegerList')
    mcols(exonReducedRanges) <- cbind(mcols(exonReducedRanges), intronReducedTable[, c('promoterId', 'MinMergedIntronRank', 'MaxMergedIntronRank', 'MaxMergedIntronTxWidth')])
    return(exonReducedRanges)
}

# Annotate all unique introns with the promoter metadata
#' @importFrom methods as
annotateIntronRanges <- function(intronRanges.unique, intronTable) {
    # Annotate unique intron ranges
    intronRanges.annotated <- intronRanges.unique
    intronRanges.annotated$MinIntronRank <- intronTable$MinIntronRank
    intronRanges.annotated$MaxIntronRank <- intronTable$MaxIntronRank
    intronRanges.annotated$MinIntronRankTxWidthMax <- intronTable$MinIntronRankTxWidthMax
    return(intronRanges.annotated)
}
