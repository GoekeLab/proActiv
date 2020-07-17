
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
  - [Quick Start](#quick-start)
  - [Creating a Promoter Annotation
    object](#creating-a-promoter-annotation-object)
  - [Estimating Promoter Activity](#estimating-promoter-activity)
  - [Limitations](#limitations)
  - [Release History](#release-history)
  - [Reference](#reference)
  - [Contributors](#contributors)

### Installation

proActiv can be installed from GitHub with:

``` r
library("devtools")
devtools::install_github("GoekeLab/proActiv")
```

### Quick Start

proActiv estimates promoter activity from RNA-Seq data. proActiv takes
as input either BAM files or junction files (TopHat2 or STAR), and a
promoter annotation object of the relevant genome. Here we demonstrate
proActiv with STAR junction files (Human genome GRCh38 Gencode v34) as
input:

``` r
library(proActiv)

files <- list.files(system.file('extdata/vignette', package = 'proActiv'), full.names = TRUE)
result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation.gencode.v34)
```

`result` is a summarizedExperiment object which can be accessed as
follows:

  - `assays(results)` returns raw/normalized promoter counts and
    absolute/relative promoter activity  
  - `metadata(results)` returns gene expression data  
  - `rowData(results)` returns a promoter to gene ID mapping and
    promoter metadata

Below we describe the steps to create a promoter annotation object and

### Creating a Promoter Annotation object

In order to quantify promoter activity, proActiv uses a set of promoters
based on genome annotations. proActiv allows the creation of a promoter
annotation object for any genome from a TxDb with the
`preparePromoterAnnotation` function. Users have the option to either
pass the file path of the GTF/GFF or TxDb to be used, or use the TxDb
object directly as input. Here, we demonstrate creating the promoter
annotation for human genome version hg 38 (Gencode v34) with both GTF
and TxDb. We use a subset of the GTF/TxDb which includes only chromosome
22 annotations.

``` r

## From GTF file path
gtf.file <- list.files(system.file('extdata/testdata/promoterAnnotation', package = 'proActiv'), pattern = 'gtf', full.names = TRUE)
promoterAnnotation.gencode.v34.chr22 <- preparePromoterAnnotation(file = gtf.file,
                                                                  species = 'Homo_sapiens')

## From TxDb object
txdb.file <- list.files(system.file('extdata/testdata/promoterAnnotation', package = 'proActiv'), pattern = 'sqlite', full.names = TRUE)
txdb <- AnnotationDbi::loadDb(txdb.file)
promoterAnnotation.gencode.v34.chr22 <- preparePromoterAnnotation(txdb = txdb, 
                                                                  species = 'Homo_sapiens')
```

proActiv provides pre-calculated promoter annotation objects for Human
genome version hg19 (Gencode v19) and hg38 (Gencode v34). The
`PromoterAnnotation` object has 3 slots:

  - `intronRanges`: Intron ranges, giving the corresponding transcripts
    of each intron  
  - `promoterIdMapping`: The id mapping between transcripts, promoter
    ids and gene ids  
  - `promoterCoordinates`: Promoter coordinates (TSS) and internal
    promoter state, along with the 3’ coordinate of the first exon

### Estimating Promoter Activity

Once promoters in the genome are identified, proActiv estimates promoter
activity at each annotated promoter from RNA-Seq data aligned with
TopHat2 or STAR. Users have the option to either pass junction files
(TopHat2 or STAR) or BAM files to `proActiv`. proActiv takes the paths
of the input files, together with the relevant promoter annotation, and
a vector describing experimental condition corresponding to each input
file, and returns a summarizedExperiment object. The returned object
summarizes promoter counts and activity, while gene expression is stored
as metadata. Row data stores promoter metadata and mean absolute
promoter activity summarized across conditions.

Below, we demonstrate running `proActiv` with input STAR junction files
and BAM files (truncated). This data is taken from the [SGNEx
project](https://github.com/GoekeLab/sg-nex-data). The reference genome
used for alignment is Gencode v34 (GRCh38). These files can can be found
at ‘extdata/vignette’:

  - extdata/vignette/SGNEx\_A549\_Illumina\_replicate1-run1.junctions.gz
  - extdata/vignette/SGNEx\_A549\_Illumina\_replicate3-run1.junctions.gz
  - extdata/vignette/SGNEx\_A549\_Illumina\_replicate5-run1.junctions.gz
  - extdata/vignette/SGNEx\_HepG2\_Illumina\_replicate2-run1.junctions.gz
  - extdata/vignette/SGNEx\_HepG2\_Illumina\_replicate4-run1.junctions.gz
  - extdata/vignette/SGNEx\_HepG2\_Illumina\_replicate5-run1.junctions.gz

<!-- end list -->

``` r
## From BAM files - genome parameter must be provided
files <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), full.names = TRUE)
result <- proActiv(files = files, 
                   promoterAnnotation = proActiv::promoterAnnotation.gencode.v34,
                   genome = 'hg38')

## From STAR junction files
files <- list.files(system.file('extdata/vignette', package = 'proActiv'), full.names = TRUE)
result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation.gencode.v34,
                   condition = rep(c('A549','HepG2'), each = 3))
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

[Demircioğlu, Deniz, et al. “A Pan-cancer Transcriptome Analysis Reveals
Pervasive Regulation through Alternative Promoters.” *Cell* 178.6
(2019):
1465-1477.](https://www.cell.com/cell/fulltext/S0092-8674\(19\)30906-7)

## Contributors

proActiv is developed and maintained by [Deniz
Demircioglu](https://github.com/dnzdmrcgl), [Joseph
Lee](https://github.com/jleechung), and [Jonathan
Göke](https://github.com/jonathangoeke).

![Stay proActiv\!](man/figures/proActiv_logoName.png)
