#' Get the reduced first exon ranges without metadata for the annotation txdb
#' provided
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#' @param numberOfCores A numeric value. The number of cores to be used for
#'   reducing first exons of each gene. Defaults to 1 (no parallelization). This
#'   parameter will be used argument to mclapply function hence require
#'   'parallel' package to be installed
#'
#' @return A GRanges object. The reduced first exon ranges for each promoter
#' @export
#'
#' @examples
#' \dontrun{
#' exonReducedRanges <- getUnannotatedReducedExonRanges(txdb,
#'                                                      species = 'Homo_sapiens',
#'                                                      numberOfCores = 1)
#' }
#'
getUnannotatedReducedExonRanges <- function(txdb, species = 'Homo_sapiens', numberOfCores = 1) {
  transcriptRanges <- getTranscriptRanges(txdb, species)
  exonRangesByTx <- getExonRangesByTx(txdb, species = species)
  transcriptRanges$TxWidth <- sapply(width(exonRangesByTx), sum)

  # Annotate first exons of each transcript with the gene/tx ids and group by gene to identify promoters
  exonRangesByTx.unlist <- unlist(exonRangesByTx)
  exonRanges.firstExon <- getFirstExonRanges(exonRangesByTx.unlist, transcriptRanges)
  exonRanges.firstExon.geneId <- as.character(transcriptRanges$GENEID)[match(exonRanges.firstExon$tx_name, transcriptRanges$TXNAME)]
  exonRanges.firstExon.byGene <- split(exonRanges.firstExon, exonRanges.firstExon.geneId)

  # Identify overlapping first exons for each gene and annotate with the promoter ids
  exonReducedRanges <- getReducedExonRanges(exonRanges.firstExon.byGene, numberOfCores)
  return(exonReducedRanges)
}

#' Prepare the id mapping between transcript ids, names, TSS ids, promoter ids
#' and gene ids
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#' @param exonReducedRanges A GRanges object. The reduced first exon ranges for
#'   each promoter
#'
#' @return A data.frame object. The id mapping between transcript ids, names,
#'   TSS ids, promoter ids and gene ids
#' @export
#'
#' @examples
#' \dontrun{
#' promoterIdMapping <- preparePromoterIdMapping(txdb,
#'                                                species = 'Homo_sapiens',
#'                                                exonReducedRanges)
#' }
#'
preparePromoterIdMapping <- function(txdb, species = 'Homo_sapiens', exonReducedRanges) {
  transcriptRanges <- getTranscriptRanges(txdb, species)
  tssCoordinates <- getTssRanges(transcriptRanges)
  exonRangesByTx <- getExonRangesByTx(txdb, species = species)
  transcriptRanges$TxWidth <- sapply(width(exonRangesByTx), sum)

  # Annotate first exons of each transcript with the gene/tx ids and group by gene to identify promoters
  exonRangesByTx.unlist <- unlist(exonRangesByTx)
  exonRanges.firstExon <- getFirstExonRanges(exonRangesByTx.unlist, transcriptRanges)
  exonRanges.firstExon.geneId <- as.character(transcriptRanges$GENEID)[match(exonRanges.firstExon$tx_name, transcriptRanges$TXNAME)]
  exonRanges.firstExon.byGene <- split(exonRanges.firstExon, exonRanges.firstExon.geneId)

  # Prepare mapping between transcripts, tss, promoters and genes
  print('Prepare mapping between transcripts, tss, promoters and genes...')
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

  return(promoterIdMapping)
}

#' Prepare the annotated intron ranges for the annotating txdb specified
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#' @param promoterIdMapping A data.frame object. The id mapping between
#'   transcript ids, names, TSS ids, promoter ids and gene ids
#'
#' @return A GRanges object. The intron ranges annotated with the promoter
#'   information.
#' @export
#'
#' @examples
#' \dontrun{
#' intronRanges.annotated <- prepareAnnotatedIntronRanges(txdb,
#'                                                        species = 'Homo_sapiens',
#'                                                        promoterIdMapping)
#' }
#'
prepareAnnotatedIntronRanges <- function(txdb, species = 'Homo_sapiens', promoterIdMapping) {
  print('Prepare annotated intron ranges...')
  transcriptRanges <- getTranscriptRanges(txdb, species)
  exonRangesByTx <- getExonRangesByTx(txdb, species)
  transcriptRanges$TxWidth <- sapply(width(exonRangesByTx), sum)

  # Retrieve intron ranges for each promoter
  intronRangesByTx <- getIntronRangesByTx(txdb, transcriptRanges, species)
  intronRangesByTx.unlist <- unlist(intronRangesByTx)
  intronRanges.unique <- getUniqueIntronRanges(intronRangesByTx.unlist)
  intronRanges.overlap <- findOverlaps(intronRangesByTx.unlist, intronRanges.unique, type = 'equal')
  intronRangesByTx.unlist$INTRONID <- subjectHits(intronRanges.overlap)

  # Annotate intron's with ranks within each transcript
  intronRankByTx <- getIntronRanks(intronRangesByTx)

  # Annotate all intron ranges with corresponding ids
  intronRangesByTx.unlist <- annotateAllIntronRanges(intronRankByTx, intronRangesByTx.unlist, transcriptRanges, promoterIdMapping)

  # Annotate unique intron ranges with corresponding ids
  intronRanges.unique <- annotateUniqueIntronRanges(intronRanges.unique, intronRangesByTx.unlist)

  intronTable <- getIntronTable(intronRanges.unique, intronRangesByTx.unlist)

  intronRanges.annotated <- annotateIntronRanges(intronRanges.unique, intronTable)

  return(intronRanges.annotated)
}

#' Prepare the reduced exon ranges for each promoter using the annotation txdb
#' specified
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#' @param promoterIdMapping A data.frame object. The id mapping between
#'   transcript ids, names, TSS ids, promoter ids and gene ids
#' @param exonReducedRanges A GRanges object. The reduced first exon ranges for
#'   each promoter without metadata
#'
#' @return A GRanges object. The reduced first exon ranges for each promoter
#'   with promoter metadata
#' @export
#'
#' @examples
#' \dontrun{
#' exonReducedRanges <- prepareAnnotatedReducedExonRanges(txdb,
#'                                                        species = 'Homo_sapiens',
#'                                                        promoterIdMapping,
#'                                                        exonReducedRanges)
#' }
#'
prepareAnnotatedReducedExonRanges <- function(txdb, species = 'Homo_sapiens', promoterIdMapping, exonReducedRanges) {
  print('Prepare the reduced exon ranges...')
  transcriptRanges <- getTranscriptRanges(txdb, species)
  exonRangesByTx <- getExonRangesByTx(txdb, species = species)
  transcriptRanges$TxWidth <- sapply(width(exonRangesByTx), sum)

  # Retrieve intron ranges for each promoter
  intronRangesByTx <- getIntronRangesByTx(txdb, transcriptRanges, species)
  intronRangesByTx.unlist <- unlist(intronRangesByTx)
  intronRanges.unique <- getUniqueIntronRanges(intronRangesByTx.unlist)
  intronRanges.overlap <- findOverlaps(intronRangesByTx.unlist, intronRanges.unique, type = 'equal')
  intronRangesByTx.unlist$INTRONID <- subjectHits(intronRanges.overlap)

  # Annotate intron's with ranks within each transcript
  intronRankByTx <- getIntronRanks(intronRangesByTx)

  # Annotate all intron ranges with corresponding ids
  intronRangesByTx.unlist <- annotateAllIntronRanges(intronRankByTx, intronRangesByTx.unlist, transcriptRanges, promoterIdMapping)

  # Annotate unique intron ranges with corresponding ids
  intronRanges.unique <- annotateUniqueIntronRanges(intronRanges.unique, intronRangesByTx.unlist)

  intronTable <- getIntronTable(intronRanges.unique, intronRangesByTx.unlist)

  intronIdByPromoter.firstIntron <- getIntronIdPerPromoter(intronRangesByTx.unlist, promoterIdMapping)

  # Retrieve first intron ranges
  intronRanges.firstIntron <- intronRangesByTx.unlist[intronRangesByTx.unlist$INTRONRANK == 1]
  intronTable.firstIntron <- intronTable[match(intronRanges.firstIntron$INTRONID, intronTable$INTRONID), ]
  intronTable.firstIntron$promoterId <- intronRanges.firstIntron$promoterId

  promoterMetadata <- getPromoterMetadata(intronTable.firstIntron)

  exonReducedRanges <- annotateReducedExonRanges(exonReducedRanges, promoterMetadata, intronIdByPromoter.firstIntron)

  return(exonReducedRanges)
}

#' Retrieve promoter coordinates (TSS) for each promoter with gene id and internal promoter state
#'
#' @param exonReducedRanges A GRanges object. The reduced first exon ranges for
#'   each promoter without metadata
#' @param promoterIdMapping A data.frame object. The id mapping between
#'   transcript ids, names, TSS ids, promoter ids and gene ids
#'
#' @return A GRanges object. Promoter coordinates (TSS) with gene id and internal promoter state
#' @export
#'
#' @examples
#' \dontrun{
#' promoterCoordinates <- preparePromoterCoordinates(exonReducedRanges,
#'                                                           promoterIdMapping)
#' }
preparePromoterCoordinates <- function(exonReducedRanges, promoterIdMapping) {
  # get TSS coordintes
  promoterCoordinates <- promoters(exonReducedRanges, downstream = 1, upstream = 0)
  promoterCoordinates.metadata <- data.frame(promoterId = promoterCoordinates$promoterId,
                                             geneId = promoterIdMapping$geneId[match(promoterCoordinates$promoterId, promoterIdMapping$promoterId)],
                                             internalPromoter = promoterCoordinates$MaxMergedIntronRank > 1,
                                             stringsAsFactors = FALSE)
  mcols(promoterCoordinates) <- promoterCoordinates.metadata
  return(promoterCoordinates)
}

#' Prepare promoter annotation data for the user specified txdb object
#'
#' @param txdb A txdb object. The txdb object of the annotation version for
#'   which promoters will be identified.
#' @param species A character object. The genus and species of the organism to
#'   be used in keepStandardChromosomes(). Supported species can be seen with
#'   names(genomeStyles()).
#' @param numberOfCores A numeric value. The number of cores to be used for
#'   reducing first exons of each gene. Defaults to 1 (no parallelization). This
#'   parameter will be used argument to mclapply function hence require
#'   'parallel' package to be installed.
#'
#' @return A PromoterAnnotation object. The reduced exon ranges, annotated
#'   intron ranges, promoter coordinates and the promoter id mapping are
#'   attributes of the promoter annotation data.
#' @export
#'
#' @examples
#' \dontrun{
#' promoterAnnotation <- preparePromoterAnnotationData(txdb,
#'                                                     species = 'Homo_sapiens',
#'                                                     numberOfCores = 1)
#' }
preparePromoterAnnotationData <- function(txdb, species = 'Homo_sapiens', numberOfCores = 1) {
  promoterAnnotationData <- PromoterAnnotation()

  # Reduce first exons to identify transcripts belonging to each promoter
  reducedExonRanges(promoterAnnotationData) <- getUnannotatedReducedExonRanges(txdb,
                                                                               species,
                                                                               numberOfCores)

  # Prepare the id mapping transcripts, TSSs, promoters and genes
  promoterIdMapping(promoterAnnotationData) <- preparePromoterIdMapping(txdb,
                                                                        species,
                                                                        reducedExonRanges(promoterAnnotationData))

  # Prepare the annotated intron ranges to be used as input for junction read counting
  annotatedIntronRanges(promoterAnnotationData) <- prepareAnnotatedIntronRanges(txdb,
                                                                                species,
                                                                                promoterIdMapping(promoterAnnotationData))

  # Annotate the reduced exons with promoter metadata
  reducedExonRanges(promoterAnnotationData) <- prepareAnnotatedReducedExonRanges(txdb,
                                                                                 species,
                                                                                 promoterIdMapping(promoterAnnotationData),
                                                                                 reducedExonRanges(promoterAnnotationData))

  # Retrieve promoter coordinates
  promoterCoordinates(promoterAnnotationData) <- preparePromoterCoordinates(reducedExonRanges(promoterAnnotationData),
                                                                            promoterIdMapping(promoterAnnotationData))

  return(promoterAnnotationData)
}
