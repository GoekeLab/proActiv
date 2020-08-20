#' Visualizes promoter activity and transcript model for a gene of interest
#'
#' @param result A SummarizedExperiment object with assays giving promoter 
#'   counts and activity with gene expression stored as column data and 
#'   promoter gene id mapping stored as row data
#' @param gene A character vector of length 1. Single gene of interest to 
#'   be plotted
#' @param txdb A TxDb object. The txdb must correspond to the genome version 
#'   used in running proActiv. Here, it is recommended to use the same txdb in 
#'   generating promoter annotations
#' @param ranges A list of GRanges. Each entry in the list should correspond to 
#'   a transcript that will be visualized, with Genomic Ranges giving the exons 
#'   corresponding to that transcript 
#' @param cex.title A numeric value. Size of axis labels. Defaults to 0.9
#' @param cex.axis A numeric value. Size of axis and axis ticks. Defaults to 0.9
#' @param cex.main A numeric value. Size of plot name. Defaults to 1
#' @param blk.width A numeric value. The width of promoters blocks in the 
#'   data track. Defaults to 500 (bases)
#' @param blk.fill A character vector of length 1. The fill colour of the 
#'   promoter blocks in the data track. Defaults to 'grey' 
#' @param blk.border A character vector of length 1. The border colour of the 
#'   promoter blocks in the data track. Defaults to 'darkgrey'
#' @param label.col A character vector of length 1. The font colour of the 
#'   promoter ID label in the annotation track. Defaults to 'black'
#' @param label.size A numeric value. The size of the promoter ID label in the 
#'   annotation track. Defaults to 0.7
#' @param arrow.width A numeric value. The width of promoter arrows in the 
#'   annotation track. This value is internally calculated based on the gene 
#'   of interest 
#' @param arrow.fill A character vector of length 1. The fill colour of the 
#'   promoter arrows in the annotation track. Defaults to 'transparent'
#' @param arrow.border A character vector of length 1. The border colour of 
#'   the promoter arrows in the annotation track. Defaults to 'grey'
#' 
#' @export
#' @return Outputs a plot of the promoters of the gene of interest across 
#'   conditions, along with a model of transcripts belonging to the gene 
#'   
#' @examples 
#'  
#' ## First, run proActiv to generate a summarizedExperiment result
#' files <- list.files(system.file('extdata/vignette/junctions', 
#'                        package = 'proActiv'), 
#'                        full.names = TRUE)
#' promoterAnnotation <- promoterAnnotation.gencode.v34.subset
#' result <- proActiv(files = files,
#'                        promoterAnnotation  = promoterAnnotation,
#'                        condition = rep(c('A549','HepG2'), each=3),
#'                        ncores = 1)
#' ## Read in pre-computed ranges
#' txdb <- AnnotationDbi::loadDb(system.file('extdata/vignette/annotations',
#'                                    'gencode.v34.annotation.rap1gap.sqlite',
#'                                    package = 'proActiv'))
#' ## Declare a gene of interest
#' gene <- 'ENSG00000076864.19'
#' ## Call plot 
#' plotPromoters(result = result, gene = gene, txdb = txdb)
#'                            
#' @importFrom Gviz plotTracks GenomeAxisTrack
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom S4Vectors complete.cases
plotPromoters <- function(result, gene, txdb, ranges,
                            cex.title = 0.9, cex.axis = 0.9, cex.main = 1,
                            blk.width = 500, blk.fill = 'grey',
                            blk.border = 'darkgrey',
                            label.col = 'black', label.size = 0.7,
                            arrow.width = NULL, arrow.fill = 'transparent', 
                            arrow.border = 'grey') {
    result.gene <- result[rowData(result)$geneId == gene, ]
    rdata <- rowData(result.gene)[complete.cases(rowData(result.gene)),]
    groups <- unique(colData(result.gene)$condition)

    if (nrow(rdata) == 0) {
        stop('Gene ID selected is either not present or has only one transcript 
            which is a single-exon transcript. proActiv does not estimate 
            promoter activity in such cases.')
    }
    print(paste0('Plotting ', gene))
    
    grtrack <- getGeneRegionTrack(rdata, gene, txdb, ranges)
    dtracklist <- getDataTrack(rdata, groups, blk.width = blk.width,
                                fill.histogram = blk.fill, 
                                col.histogram = blk.border)
    atrack <- getAnnotationTrack(rdata, grtrack, arrow.width = arrow.width,
                                fontcolor.feature = label.col,
                                cex.feature = label.size,
                                fill = arrow.fill, col = arrow.border)
    gtrack <- GenomeAxisTrack()
    
    print('Creating Plot...')
    plotTracks(c(grtrack, dtracklist, atrack, gtrack), type='histo', 
                main = gene, col.axis = 'black', col.title = 'black', 
                background.title = 'transparent', cex.title = cex.title, 
                cex.axis = cex.axis, cex.main = cex.main)
    }

# Helper function to get data track
#' @importFrom GenomicRanges GRanges values
#' @importFrom Gviz DataTrack
getDataTrack <- function(rdata, groups, blk.width,
                            fill.histogram, col.histogram) {
    print('Creating Data Track...')
    ## Set dummy width for visualization
    if (is.null(blk.width)) {
        blk.width <- blk.width
    }
    blk.offset <- 0
    if (unique(rdata$strand) == '-') {
        blk.offset <- -blk.width
    } 
    ## GRanges summarizing promoter activity and coordinates
    gr <- GRanges(seqnames = rdata$seqnames, 
                    strand = rdata$strand,
                    ranges = IRanges(start = rdata$start + blk.offset, 
                                        width = blk.width)) 
    values(gr) <- rdata[,grep('mean', colnames(rdata))]
    names(mcols(gr)) <- gsub('.mean', '', names(mcols(gr)))

    ## Set max y limit
    maxy <- max(unlist(mcols(gr))) + 0.5
    dtracklist <- vector("list", length(groups))
    for (i in seq_len(length(groups))) {
        dtracklist[[i]] <- DataTrack(gr[,groups[i]], 
                                    name = groups[i], ylim = c(0, maxy),
                                    fill.histogram = fill.histogram,
                                    col.histogram = col.histogram)
    }
    return(dtracklist)
    }

# Helper function to get gene region track
#' @importFrom GenomicFeatures exonsBy
#' @importFrom Gviz GeneRegionTrack
#' @importFrom GenomicRanges GRangesList
getGeneRegionTrack <- function(rdata, gene, txdb, ranges) {
    print('Creating Gene Region Track...')
    txs.gene <- rdata$txId[[1]]
    if ( missing(txdb) & missing(ranges)) {
        stop('Either txdb or ranges must be provided')
    } else if (!missing(txdb) & !missing(ranges)) {
        stop('Both `txdb` and `ranges` arguments are provided: 
                Please specify only one')
    } else if (missing(ranges)) {
        exons.gene <- suppressWarnings(exonsBy(txdb, by = 'tx', 
                                               use.names = TRUE)[txs.gene])
    } else {
        exons.gene <- ranges
    }
    gene.model <- as(GRangesList(exons.gene), 'data.frame')
    gene.model <- gene.model[,c('seqnames','start','end'
                                ,'strand','group_name')]
    colnames(gene.model)[1] <- 'chromosome'
    colnames(gene.model)[5] <- 'transcript'
    
    grtrack <- GeneRegionTrack(gene.model,
                            transcriptAnnotation = 'transcript',
                            name = as.character(gene.model$chromosome[1]))
    return(grtrack)
    }

# Helper function to get annotation track
#' @importFrom Gviz AnnotationTrack
getAnnotationTrack <- function(rdata, grtrack, arrow.width,
                                fontcolor.feature, cex.feature, 
                                fill, col) {
    print('Creating Annotation Track...')
    ## Adjust arrow start depending on stand - Gviz quirks
    if (is.null(arrow.width)) {
        min.coord <- min(c(start(grtrack), end(grtrack)))
        max.coord <- max(c(start(grtrack), end(grtrack)))
        arrow.width <- diff(c(min.coord, max.coord))/sqrt(nrow(rdata)+1)
    }
    arrow.offset <- 0
    if (unique(rdata$strand) == '-') {
        arrow.offset <- -arrow.width
    } 
    atrack <- AnnotationTrack(start = rdata$start + arrow.offset, 
                            width = arrow.width,
                            chromosome = rdata$seqnames, 
                            strand = rdata$strand,
                            group = paste0('prmtr.', rdata$promoterId),
                            featureAnnotation = 'group', name = '',
                            fontcolor.feature = fontcolor.feature,
                            cex.feature = cex.feature,
                            fill = fill, col = col)
    return(atrack)
    }
