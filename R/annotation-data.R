#' Prepares promoter annotation from a gtf or txdb 
#'
#' @param txdb A txdb object. The txdb of the annotation version for which
#'   promoters will be identified. Either `txdb` or `file` argument must be 
#'   specified, but not both.
#' @param file A character object. The file path to a gtf/gff or txdb of the 
#'   annotation version for which promoters will be identified. Either `txdb` 
#'   or `file` argument must be specified, but not both.
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#'
#' @return A PromoterAnnotation object. The annotated
#'   intron ranges, promoter coordinates and the promoter id mapping are
#'   attributes of the promoter annotation data.
#' @export
#'
#' @examples
#' 
#' txdbPath <- system.file('extdata/vignette/annotations/',
#'                             'gencode.v34.annotation.subset.sqlite', 
#'                              package = 'proActiv')
#' txdb <- AnnotationDbi::loadDb(txdbPath)
#' promoterAnnotation <- preparePromoterAnnotation(txdb = txdb,
#'                                                 species = 'Homo_sapiens')
#' 
preparePromoterAnnotation <- function(txdb, file, species) {
    txdb <- parseTxdb(txdb, file, species)
    promoterAnnotation <- PromoterAnnotation()
    
    ### Reduce first exons to identify transcripts belonging to each promoter 
    print('Extract exons by transcripts...')
    transcriptRanges <- getTranscriptRanges(txdb, species) 
    exonRangesByTx <- getExonRangesByTx(txdb, species = species) 
    transcriptRanges$TxWidth <- vapply(width(exonRangesByTx), 
                                        FUN = sum, FUN.VALUE = numeric(1)) 
    
    ## Annotate first exons of each transcript 
    exonRangesByTx.unlist <- unlist(exonRangesByTx)
    exonRanges.firstExon <- getFirstExonRanges(exonRangesByTx.unlist)
    exonRanges.firstExon.geneId <- as.character(
        transcriptRanges$GENEID)[match(exonRanges.firstExon$tx_name, 
                                        transcriptRanges$TXNAME)]
    
    ## Identify overlapping first exons by gene and annotate with promoter ids
    print('Identify overlapping first exons for each gene...')
    exonReducedRanges <- getReducedExonRanges(exonRanges.firstExon, 
                                                exonRanges.firstExon.geneId)
    
    ### Prepare the id mapping transcripts, TSSs, promoters and genes ###
    print('Prepare mapping between transcripts, tss, promoters and genes...')
    tssCoordinates <- getTssRanges(transcriptRanges)
    revmapTMP <- as.vector(exonReducedRanges$revmap)
    txName.all <- exonRanges.firstExon$tx_name[unlist(revmapTMP)]
    promoterId.all <- rep(exonReducedRanges$promoterId, vapply(revmapTMP, 
                                                        FUN = length, 
                                                        FUN.VALUE = numeric(1)))
    txId.all <- transcriptRanges$TXID[match(txName.all, transcriptRanges$TXNAME)]
    tssId.all <- tssCoordinates$tssId[match(txName.all, tssCoordinates$TXNAME)]
    geneId.all <- as.character(transcriptRanges$GENEID)[match(txName.all, transcriptRanges$TXNAME)]
    promoterIdMapping <- data.frame(txId = txId.all,
                                    txName = txName.all,
                                    tssId = tssId.all,
                                    promoterId = promoterId.all,
                                    geneId = geneId.all,
                                    stringsAsFactors = FALSE)
    promoterIdMapping <- promoterIdMapping[match(seq_len(nrow(promoterIdMapping)), promoterIdMapping$txId), ]
    rownames(promoterIdMapping) <- NULL
    colnames(promoterIdMapping) <- c('transcriptId', 'transcriptName', 
                                        'tssId', 'promoterId', 'geneId')
    
    ### Prepare the annotated intron ranges to be used as input for junction read counting ###
    print('Prepare annotated intron ranges...')
    intronRangesByTx <- getIntronRangesByTx(txdb, transcriptRanges, species) 
    intronRangesByTx.unlist <- unlist(intronRangesByTx)
    intronRanges.unique <- getUniqueIntronRanges(intronRangesByTx.unlist) 
    
    ## Annotate full length intron ranges with intron id
    intronRanges.overlap <- findOverlaps(intronRangesByTx.unlist, 
                                        intronRanges.unique, type = 'equal')
    intronRangesByTx.unlist$INTRONID <- subjectHits(intronRanges.overlap)
    
    ## Extract ranks of introns for each transcript
    intronRankByTx <- getIntronRanks(intronRangesByTx)
    ## Annotate all intron ranges with corresponding ids
    intronRangesByTx.unlist <- annotateAllIntronRanges(intronRankByTx, 
                                                        intronRangesByTx.unlist, 
                                                        transcriptRanges, 
                                                        promoterIdMapping)
    
    ## Annotate unique intron ranges with corresponding ids
    intronRanges.unique <- annotateUniqueIntronRanges(intronRanges.unique, intronRangesByTx.unlist)
    intronTable <- getIntronTable(intronRanges.unique, intronRangesByTx.unlist)
    intronRanges.annotated <- annotateIntronRanges(intronRanges.unique, intronTable)
    
    ### Annotate the reduced exons with promoter metadata
    print('Annotating reduced exon ranges...')
    ## Identify first intron for each promoter for each transcript 
    intronIdByPromoter.firstIntron <- getIntronIdPerPromoter(intronRangesByTx.unlist, promoterIdMapping)
    
    ## Retrieve first intron ranges
    intronRanges.firstIntron <- intronRangesByTx.unlist[intronRangesByTx.unlist$INTRONRANK == 1]
    intronTable.firstIntron <- intronTable[match(intronRanges.firstIntron$INTRONID, intronTable$INTRONID), ]
    intronTable.firstIntron$promoterId <- intronRanges.firstIntron$promoterId
    
    promoterMetadata <- getPromoterMetadata(intronTable.firstIntron)
    exonReducedRanges <- annotateReducedExonRanges(exonReducedRanges, promoterMetadata, intronIdByPromoter.firstIntron)
    
    ### Retrieve promoter coordinates
    print('Prepare promoter coordinates and first exon ranges...')
    promoterCoordinates <- promoters(exonReducedRanges, downstream = 1, upstream = 0)
    exonEnd <- resize(exonReducedRanges, width = 1, fix = 'end')
    promoterCoordinates.metadata <- data.frame(promoterId = promoterCoordinates$promoterId,
                                                geneId = promoterIdMapping$geneId[match(promoterCoordinates$promoterId, promoterIdMapping$promoterId)],
                                                internalPromoter = promoterCoordinates$MaxMergedIntronRank > 1,
                                                firstExonEnd = end(exonEnd),
                                                stringsAsFactors = FALSE)
    mcols(promoterCoordinates) <- promoterCoordinates.metadata
    promoterCoordinates$intronId <- exonReducedRanges$intronId
    
    ### Build Promoter Annotation
    intronRanges(promoterAnnotation) <- intronRanges.annotated[,c('INTRONID', 'TXNAME')]
    promoterIdMapping(promoterAnnotation) <- promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')]
    promoterCoordinates(promoterAnnotation) <- promoterCoordinates[,c('promoterId', 'internalPromoter', 'firstExonEnd', 'intronId')]
    return(promoterAnnotation)
} 
