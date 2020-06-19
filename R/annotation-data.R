#' Prepare promoter annotation data for the user specified txdb object
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified.
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#'
#' @return A PromoterAnnotation object. The reduced exon ranges, annotated
#'   intron ranges, promoter coordinates and the promoter id mapping are
#'   attributes of the promoter annotation data.
#' @export
#'
#' @examples
#' \dontrun{
#' promoterAnnotation <- preparePromoterAnnotation(txdb,
#'                                                 species = 'Homo_sapiens')
#' }
preparePromoterAnnotation <- function(txdb, species) {
  
  ## Instantiate Promoter Annotation
  promoterAnnotation <- PromoterAnnotation()
  
  ############################################################################# 
  ### Reduce first exons to identify transcripts belonging to each promoter ###
  ############################################################################# 
  
  ## Extract Transcript Ranges
  print('Extract exons by transcripts...')
  transcriptRanges <- getTranscriptRanges(txdb, species) 
  exonRangesByTx <- getExonRangesByTx(txdb, species = species) 
  transcriptRanges$TxWidth <- sapply(width(exonRangesByTx), sum) 
  
  ## Annotate first exons of each transcript 
  exonRangesByTx.unlist <- unlist(exonRangesByTx)
  exonRanges.firstExon <- getFirstExonRanges(exonRangesByTx.unlist)
  exonRanges.firstExon.geneId <- as.character(transcriptRanges$GENEID)[match(exonRanges.firstExon$tx_name, transcriptRanges$TXNAME)]
  
  ## Identify overlapping first exons for each gene and annotate with promoter ids
  print('Identify overlapping first exons for each gene...')
  exonReducedRanges <- getReducedExonRanges(exonRanges.firstExon, exonRanges.firstExon.geneId)
  
  ##################################################################### 
  ### Prepare the id mapping transcripts, TSSs, promoters and genes ###
  ##################################################################### 
  
  print('Prepare mapping between transcripts, tss, promoters and genes...')
  tssCoordinates <- getTssRanges(transcriptRanges)
  revmapTMP <- as.vector(exonReducedRanges$revmap)
  txName.all <- exonRanges.firstExon$tx_name[unlist(revmapTMP)]
  promoterId.all <- rep(exonReducedRanges$promoterId, sapply(revmapTMP, length))
  txId.all <- transcriptRanges$TXID[match(txName.all, transcriptRanges$TXNAME)]
  tssId.all <- tssCoordinates$tssId[match(txName.all, tssCoordinates$TXNAME)]
  geneId.all <- as.character(transcriptRanges$GENEID)[match(txName.all, transcriptRanges$TXNAME)]
  promoterIdMapping <- data.frame(txId = txId.all,
                                  txName = txName.all,
                                  tssId = tssId.all,
                                  promoterId = promoterId.all,
                                  geneId = geneId.all,
                                  stringsAsFactors = FALSE)
  promoterIdMapping <- promoterIdMapping[match(1:nrow(promoterIdMapping), promoterIdMapping$txId), ]
  rownames(promoterIdMapping) <- NULL
  colnames(promoterIdMapping) <- c('transcriptId', 'transcriptName', 'tssId', 'promoterId', 'geneId')
  
  ########################################################################################## 
  ### Prepare the annotated intron ranges to be used as input for junction read counting ###
  ##########################################################################################
  
  print('Prepare annotated intron ranges...')
  intronRangesByTx <- getIntronRangesByTx(txdb, transcriptRanges, species) 
  intronRangesByTx.unlist <- unlist(intronRangesByTx)
  intronRanges.unique <- getUniqueIntronRanges(intronRangesByTx.unlist) 
  
  ## Annotate full length intron ranges with intron id
  intronRanges.overlap <- findOverlaps(intronRangesByTx.unlist, intronRanges.unique, type = 'equal')
  intronRangesByTx.unlist$INTRONID <- subjectHits(intronRanges.overlap)
  
  ## Extract ranks of introns for each transcript
  intronRankByTx <- getIntronRanks(intronRangesByTx)
  
  ## Annotate all intron ranges with corresponding ids
  intronRangesByTx.unlist <- annotateAllIntronRanges(intronRankByTx, intronRangesByTx.unlist, transcriptRanges, promoterIdMapping) #OK
  
  ## Annotate unique intron ranges with corresponding ids
  intronRanges.unique <- annotateUniqueIntronRanges(intronRanges.unique, intronRangesByTx.unlist)
  intronTable <- getIntronTable(intronRanges.unique, intronRangesByTx.unlist)
  intronRanges.annotated <- annotateIntronRanges(intronRanges.unique, intronTable)
  
  #########################################################
  ### Annotate the reduced exons with promoter metadata ###
  #########################################################
  
  print('Prepare reduced exon ranges...')
  ## Identify first intron for each promoter for each transcript 
  intronIdByPromoter.firstIntron <- getIntronIdPerPromoter(intronRangesByTx.unlist, promoterIdMapping) # this replaces revmap
  
  ## Retrieve first intron ranges
  intronRanges.firstIntron <- intronRangesByTx.unlist[intronRangesByTx.unlist$INTRONRANK == 1]
  intronTable.firstIntron <- intronTable[match(intronRanges.firstIntron$INTRONID, intronTable$INTRONID), ]
  intronTable.firstIntron$promoterId <- intronRanges.firstIntron$promoterId
  
  promoterMetadata <- getPromoterMetadata(intronTable.firstIntron)
  
  exonReducedRanges <- annotateReducedExonRanges(exonReducedRanges, promoterMetadata, intronIdByPromoter.firstIntron)
  
  #####################################
  ### Retrieve promoter coordinates ###
  #####################################
  
  print('Prepare promoter coordinates...')
  promoterCoordinates <- promoters(exonReducedRanges, downstream = 1, upstream = 0)
  promoterCoordinates.metadata <- data.frame(promoterId = promoterCoordinates$promoterId,
                                             geneId = promoterIdMapping$geneId[match(promoterCoordinates$promoterId, promoterIdMapping$promoterId)],
                                             internalPromoter = promoterCoordinates$MaxMergedIntronRank > 1,
                                             stringsAsFactors = FALSE)
  mcols(promoterCoordinates) <- promoterCoordinates.metadata
  
  #################################
  ### Build Promoter Annotation ###
  #################################
  
  reducedExonRanges(promoterAnnotation) <- exonReducedRanges
  annotatedIntronRanges(promoterAnnotation) <- intronRanges.annotated
  promoterIdMapping(promoterAnnotation) <- promoterIdMapping
  promoterCoordinates(promoterAnnotation) <- promoterCoordinates
  return(promoterAnnotation)
} 

