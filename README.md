
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

proActiv estimates promoter activity from RNA-Seq data. Promoter
activity is defined as the total amount of transcription initiated at
each promoter. proActiv takes as input either BAM files or junction
files (TopHat2 or STAR), and a promoter annotation object of the
relevant genome. An optional argument `condition` can be supplied,
describing the condition corresponding to each input file. Here we
demonstrate proActiv with STAR junction files (Human genome GRCh38
GENCODE v34) as input:

``` r
library(proActiv)

## List of STAR junction files as input
files <- list.files(system.file('extdata/vignette/junctions', 
                                package = 'proActiv'), full.names = TRUE)
## Vector describing experimental condition
condition <- rep(c('A549','HepG2'), each=3)
## Promoter annotation for human genome GENCODE v34
promoterAnnotation <- promoterAnnotation.gencode.v34

result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation,
                   condition = condition)
```

`result` is a summarizedExperiment object which can be accessed as
follows:

  - `assays(results)` returns raw/normalized promoter counts and
    absolute/relative promoter activity  
  - `metadata(results)` returns gene expression data  
  - `rowData(results)` returns promoter metadata and summarized absolute
    promoter activity by conditions

### Creating a Promoter Annotation object

In order to quantify promoter activity, proActiv uses a set of promoters
based on genome annotations. proActiv allows the creation of a promoter
annotation object for any genome from a TxDb object or from a GTF file
with the `preparePromoterAnnotation` function. Users have the option to
either pass the file path of the GTF/GFF or TxDb to be used, or use the
TxDb object directly as input. proActiv includes pre-calculated promoter
annotations for human and mouse genomes:

  - GENCODE Release 19 / GRCh37 / hg19 :
    `promoterAnnotation.gencode.v19`
  - GENCODE Release 34 / GRCh38 / hg38 :
    `promoterAnnotation.gencode.v34`
  - GENCODE Release M1 / NCBIM37 / mm9 :
    `promoterAnnotation.gencode.vM1`
  - GENCODE Release M25 / GRCm38.p6 / mm10 :
    `promoterAnnotation.gencode.vM25`

GTF files can be downloaded from the
[GENCODE](https://www.gencodegenes.org) page.

Here, we demonstrate creating the promoter annotation for the Human
genome (GENCODE v34) with both GTF and TxDb. To keep the run-time small,
we use a subset of the GTF/TxDb which includes only chromosome 22
annotations.

``` r
## From GTF file path
gtf.file <- system.file('extdata/vignette/annotation/gencode.v34.annotation.chr22.gtf.gz', 
                        package = 'proActiv')
promoterAnnotation.gencode.v34.chr22 <- preparePromoterAnnotation(file = gtf.file,
                                                                  species = 'Homo_sapiens')
## From TxDb object
txdb.file <- system.file('extdata/vignette/annotation/gencode.v34.annotation.chr22.sqlite', 
                         package = 'proActiv')
txdb <- loadDb(txdb.file)
promoterAnnotation.gencode.v34.chr22 <- preparePromoterAnnotation(txdb = txdb, 
                                                                  species = 'Homo_sapiens')
```

The `PromoterAnnotation` object has 3 slots:

  - `intronRanges` : Intron ranges, giving the corresponding transcripts
    of each intron
  - `promoterIdMapping` : An ID mapping between transcripts, promoter
    IDs and gene IDs  
  - `promoterCoordinates` : Promoter coordinates (TSS) and internal
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
used for alignment is GENCODE v34 (GRCh38). These files can can be found
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
                   promoterAnnotation = promoterAnnotation.gencode.v34,
                   genome = 'hg38')

## From STAR junction files
files <- list.files(system.file('extdata/vignette', package = 'proActiv'), full.names = TRUE)
result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation.gencode.v34,
                   condition = rep(c('A549','HepG2'), each = 3))
```

## Release History

**Release 0.99.0**

Release date: 4th August 2020

Changes in version 0.99.0:

  - Workflow: The wrapper function `proActiv` performs all steps to
    estimate promoter activity and calculates promoter metadata. A
    `condition` argument can be supplied for `proActiv` to summarize
    promoter counts and activity across conditions. These results are
    returned as a SummarizedExperiment object.

  - BAM file usage: In addition to junction files, `proActiv` now allows
    BAM files as input. However, users should note that this function is
    not fully optimized and may have long run-time.

  - Promoter annotation: Improved efficiency in generating promoter
    annotations without the need for parallelization with the
    `preparePromoterAnnotations` function. Promoter annotation objects
    for human (hg19/hg38) and mouse (mm9/mm10) genomes are now
    pre-calculated and available to the user. The promoter annotation
    object is also trimmed to preserve essential information for running
    `proActiv`, in order to comply with Bioconductor guidelines
    concerning package size.

  - Plotting promoter activity: The plotting function `plotPromoters`
    visualizes promoter activity across conditions. It accepts the
    SummarizedExperiment object returned by `proActiv` along with a gene
    of interest and gene annotations as arguments. This allows users to
    visualize promoter activity and identify instances of alternative
    promoter usage.

  - Vignette: proActiv now comes with a vignette, documenting a complete
    step-by-step workflow in identifying active and alternative promoter
    usage. This includes guidance on running `proActiv`, creating
    promoter annotations and identifying alternative promoter usage.
    Various visualizations of promoter activity are also offered.

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
