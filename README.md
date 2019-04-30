
<!-- README.md is generated from README.Rmd. Please edit that file -->
proActiv
--------

<!-- badges: start -->
<!-- badges: end -->
proActiv is an R package that estimates promoter activity from RNA-Seq data. proActiv uses aligned reads and genome annotations as input, and provides absolute and relative promoter activity as output. The package can be used to identify active promoters and alternative promoters, the details of the method are described at <https://doi.org/10.1101/176487>.

### Installation

proActiv can be installed from [GitHub](https://github.com/) with:

``` r
library("devtools")
devtools::install_github("GoekeLab/proActiv")
```

### Annotation and Example Data

The pre-calculated promoter annotation objects can be downloaded here for hg19 (Grch37): [preprocessedAnnotation](https://drive.google.com/drive/folders/1dtuP2QIKBTIQd8HecUmDZz2fCI_VaLMl?usp=sharing)

Example junction files as produced by TopHat2 and STAR can be downloaded here: [inputFiles](https://drive.google.com/drive/folders/1R8sI97h1ZTdyxbQxG4latR9xN9FF2tq8?usp=sharing)

### Estimate Promoter Activity (TopHat2 alignment)

This is a basic example to estimate promoter activity from a set of RNA-Seq data which was aligned with TopHat2. proActiv will use the junction file from the TopHat2 alignment (see below for an example with STAR-aligned reads), and a set of annotation objects that describe the associations of promoters, transcripts, and genes, to calculate promoter activity.

``` r
library(proActiv)

# The number of cores to be used for parallel execution (mc.cores argument for parallel::mclappy), optional
numberOfCores <- 1

# Loads the promoter annotations for the human genome (hg19):
# exonReducedRanges, promoterIdMapping, intronRanges.annotated and promoterCoordinates
load('preprocessedAnnotation.RData')

### TopHat2 Junction Files Example 

# The paths and labels for samples
tophatJunctionFiles <- list.files('tophat/', full.names = TRUE)
tophatJunctionFileLabels <- paste0('s', 1:length(tophatJunctionFiles), '-tophat')

# Count the total number of junction reads for each promoter
promoterCounts.tophat <- calculatePromoterReadCounts(exonReducedRanges, 
                                                      intronRanges.annotated, 
                                                      junctionFilePaths = tophatJunctionFiles, 
                                                      junctionFileLabels =  tophatJunctionFileLabels, 
                                                      junctionType = 'tophat', 
                                                      numberOfCores)

# Normalize promoter read counts by DESeq2 (optional)
normalizedPromoterCounts.tophat <- normalizePromoterReadCounts(promoterCounts.tophat)

# Calculate absolute promoter activity
absolutePromoterActivity.tophat <- getAbsolutePromoterActivity(normalizedPromoterCounts.tophat, 
                                                               promoterIdMapping, 
                                                               log2 = TRUE, 
                                                               pseudocount = 1)
# Calculate gene expression
geneExpression.tophat <- getGeneExpression(absolutePromoterActivity.tophat)
# Calculate relative promoter activity
relativePromoterActivity.tophat <- getRelativePromoterActivity(absolutePromoterActivity.tophat, 
                                                               geneExpression.tophat)
```

### Estimate Promoter Activity (STAR alignment)

``` r
library(proActiv)

# The number of cores to be used for parallel execution (mc.cores argument for parallel::mclappy), optional
numberOfCores <- 1

# Loads the promoter annotations for the human genome (hg19):
# exonReducedRanges, promoterIdMapping, intronRanges.annotated and promoterCoordinates
load('preprocessedAnnotation.RData')

### STAR Junction Files Example 

# The paths and labels for samples
starJunctionFiles <- list.files('star/', full.names = TRUE)
starJunctionFileLabels <- paste0('s', 1:length(starJunctionFiles), '-star')

# Count the total number of junction reads for each promoter
promoterCounts.star <- calculatePromoterReadCounts(exonReducedRanges, 
                                                      intronRanges.annotated, 
                                                      junctionFilePaths = starJunctionFiles, 
                                                      junctionFileLabels =  starJunctionFileLabels, 
                                                      junctionType = 'star', 
                                                      numberOfCores)

# Normalize promoter read counts by DESeq2 (optional)
normalizedPromoterCounts.star <- normalizePromoterReadCounts(promoterCounts.star)

# Calculate absolute promoter activity
absolutePromoterActivity.star <- getAbsolutePromoterActivity(normalizedPromoterCounts.star, 
                                                             promoterIdMapping, 
                                                             log2 = TRUE, 
                                                             pseudocount = 1)
# Calculate gene expression
geneExpression.star <- getGeneExpression(absolutePromoterActivity.star)
# Calculate relative promoter activity
relativePromoterActivity.star <- getRelativePromoterActivity(absolutePromoterActivity.star, 
                                                             geneExpression.tophat)
```

### Creating your own promoter annotations

proActiv provides functions to create promoter annotation objects for any genome. Here we describe how the annotation can be created using a TxDb object (please see the TxDb documentation for how to create annotations from a GTF file).

A TxDb object for the human genome version hg19 (Grch37) can be downloaded here: [inputFiles](https://drive.google.com/drive/folders/1R8sI97h1ZTdyxbQxG4latR9xN9FF2tq8?usp=sharing)

``` r
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(dplyr)
library(proActiv)

# Load the txdb object for your annotation of choice (Gencode v19 used here)
txdb <- loadDb('./inputFiles/annotation/gencode.v19.annotation.sqlite')

# The species argument to be used for GenomeInfoDb::keepStandardChromosomes
species <- 'Homo_sapiens'
# The number of cores to be used for parallel execution (mc.cores argument for parallel::mclappy), optional
numberOfCores <- 1

### Annotation data preparation
### Needs to be executed once per annotation. Results can be saved and loaded later for reuse 

# Reduce first exons to identify transcripts belonging to each promoter
exonReducedRanges <- getUnannotatedReducedExonRanges(txdb, 
                                                     species,
                                                     numberOfCores)

# Prepare the id mapping transcripts, TSSs, promoters and genes
promoterIdMapping <- preparePromoterIdMapping(txdb, 
                                              species,
                                              exonReducedRanges)

# Prepare the annotated intron ranges to be used as input for junction read counting
intronRanges.annotated <- prepareAnnotatedIntronRanges(txdb, 
                                                        species, 
                                                        promoterIdMapping)

# Annotate the reduced exons with promoter metadata
exonReducedRanges <- prepareAnnotatedReducedExonRanges(txdb, 
                                                       species, 
                                                       promoterIdMapping, 
                                                       exonReducedRanges)

# Retrieve promoter coordinates 
promoterCoordinates <- preparePromoterCoordinates(exonReducedRanges,
                                                    promoterIdMapping)
```

Limitations
-----------

proActiv will not provide promoter activity estimates for promoters which are not uniquely identifiable from splice junctions (single exon transcripts, promoters which overlap with internal exons).

Citing proActiv
---------------

If you use proActiv, please cite: Demircioğlu, Deniz, et al. "A Pan-Cancer Transcriptome Analysis Reveals Pervasive Regulation through Tumor-Associated Alternative Promoters." bioRxiv (2018): 176487.

Contributors
------------

ProActiv is developed and maintained by Deniz Demircioglu and Jonathan Göke.
