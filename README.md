
<!-- README.md is generated from README.Rmd. Please edit that file -->

![Stay proActiv\!](man/figures/proActiv_design.png)

![Stay proActiv\!](man/figures/proActiv_name.png)

## proActiv: Estimation of Promoter Activity from RNA-Seq data

<!-- badges: start -->

<!-- badges: end -->

proActiv is an R package that estimates promoter activity from RNA-Seq
data. proActiv uses aligned reads and genome annotations as input, and
provides absolute and relative promoter activity as output. The package
can be used to identify active promoters and alternative promoters, the
details of the method are described in [Demircioglu et al](#reference).

Additional data on differential promoters in tissues and cancers from
TCGA, ICGC, GTEx, and PCAWG can be downloaded here:
<https://jglab.org/data-and-software/>

### Content

  - [Installation](#installation)
  - [Estimate Promoter Activity](#estimate-promoter-activity)
  - [Annotation and Example Data](#annotation-and-example-data)
  - [Creating your own promoter
    annotations](#creating-your-own-promoter-annotations)
  - [Limitations](#limitations)
  - [Release History](#release-history)
  - [Citing proActiv](#reference)
  - [Contributors](#contributors)

### Installation

proActiv can be installed from GitHub with:

``` r
library("devtools")
devtools::install_github("GoekeLab/proActiv")
```

### Estimate Promoter Activity

This is a basic example to estimate promoter activity from a set of
RNA-Seq data which was aligned with TopHat2 (or STAR). proActiv will use
the junction file from the TopHat2 (STAR) alignment, and a set of
annotation objects that describe the associations of promoters,
transcripts, and genes, to calculate promoter activity.

``` r
library(proActiv)

# Preprocessed annotations are available as part of the R package for the human genome (hg19):
# proActiv::promoterAnnotationData.gencode.v19

# The paths and labels for samples
junctionFiles <- list.files(system.file('extdata/tophat2', package = 'proActiv'), full.names = TRUE)

# for STAR alignment
# junctionFiles <- list.files(system.file('extdata/star', package = 'proActiv'), full.names = TRUE)

junctionFileLabels <- paste0('s', 1:length(junctionFiles))

# Count the total number of junction reads for each promoter
promoterCounts <- calculatePromoterReadCounts(proActiv::promoterAnnotationData.gencode.v19,
                                                      junctionFilePaths = junctionFiles,
                                                      junctionFileLabels =  junctionFileLabels,
                                                      junctionType = 'tophat')  # use junctionType = 'star' for STAR aligned reads

# Normalize promoter read counts by DESeq2 (optional)
normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts)

# Calculate absolute promoter activity
absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts,
                                                               proActiv::promoterAnnotationData.gencode.v19)
# Calculate gene expression
geneExpression <- getGeneExpression(absolutePromoterActivity)

# Calculate relative promoter activity
relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity,
                                                               geneExpression)
```

### Annotation and Example Data

Pre-calculated promoter annotation data for Gencode v19 (GRCh37) is
available as part of the proActiv package. The PromoterAnnotation object
has 4 slots:

  - reducedExonRanges : The reduced first exon ranges for each promoter
    with promoter metadata for Gencode v19
  - promoterIdMapping : The id mapping between transcript ids, names,
    TSS ids, promoter ids and gene ids for Gencode v19
  - annotatedIntronRanges : The intron ranges annotated with the
    promoter information for Gencode v19
  - promoterCoordinates : Promoter coordinates (TSS) with gene id and
    internal promoter state for Gencode v19

Example junction files as produced by TopHat2 and STAR are available as
external data. The reference genome used for alignment is Gencode v19
(GRCh37). The TopHat2 and STAR example files (5 files each) can be found
at ‘extdata/tophat2’ and ‘extdata/star’ folders respectively.

Example TopHat2 files:

  - extdata/tophat2/sample1.bed
  - extdata/tophat2/sample2.bed
  - extdata/tophat2/sample3.bed
  - extdata/tophat2/sample4.bed
  - extdata/tophat2/sample5.bed

Example STAR files:

  - extdata/tophat2/sample1.junctions
  - extdata/tophat2/sample2.junctions
  - extdata/tophat2/sample3.junctions
  - extdata/tophat2/sample4.junctions
  - extdata/tophat2/sample5.junctions

### Creating your own promoter annotations

proActiv provides functions to create promoter annotation objects for
any genome. Here we describe how the annotation can be created using a
TxDb object (please see the TxDb documentation for how to create
annotations from a GTF file).

A TxDb object for the human genome version hg19 (Grch37) can be
downloaded here:
[inputFiles](http://s3.ap-southeast-1.amazonaws.com/all-public-data.store.genome.sg/DemirciogluEtAl2019/annotations/gencode.v19.annotation.sqlite)

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
promoterAnnotationData <- preparePromoterAnnotationData(txdb, species = species, numberOfCores = numberOfCores)

# Retrieve the id mapping between transcripts, TSSs, promoters and genes
head(promoterIdMapping(promoterAnnotationData))

# Retrieve promoter coordinates
head(promoterCoordinates(promoterAnnotationData))
```

## Release History

**Initial Release 0.1.0**

Release date: 19th May 2020

This release corresponds to the proActiv version used by [Demircioglu et
al.](#reference)

## Limitations

proActiv will not provide promoter activity estimates for promoters
which are not uniquely identifiable from splice junctions (single exon
transcripts, promoters which overlap with internal exons).

## Reference

If you use proActiv, please cite:

[Demircioğlu, Deniz, et al. “A Pan-cancer Transcriptome Analysis Reveals
Pervasive Regulation through Alternative Promoters.” *Cell* 178.6
(2019):
1465-1477.](https://www.cell.com/cell/fulltext/S0092-8674\(19\)30906-7)

## Contributors

proActiv is developed and maintained by [Deniz
Demircioglu](https://github.com/dnzdmrcgl), [Joseph
Lee](https://github.com/jleechung), and [Jonathan
Göke](https://github.com/jonathangoeke).

![Stay proActiv\!](man/figures/proActiv_logoName.png)
