
<!-- README.md is generated from README.Rmd. Please edit that file -->

![Stay proActiv\!](man/figures/proActiv_design.png)

![Stay proActiv\!](man/figures/proActiv_name.png)

## proActiv: Estimation of Promoter Activity from RNA-Seq data

<!-- badges: start -->

[![GitHub release (latest by
date)](https://img.shields.io/github/v/release/GoekeLab/proActiv)](https://github.com/GoekeLab/proActiv/releases/)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-Yes-brightgreen)](https://github.com/GoekeLab/proActiv/graphs/contributors)
[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
<!-- badges: end -->

proActiv is an R package that estimates promoter activity from RNA-Seq
data. proActiv uses aligned reads and genome annotations as input, and
provides absolute and relative promoter activity as output. The package
can be used to identify active promoters and alternative promoters.
Details of the method are described in [Demircioglu et al](#reference).

HTML documentation of proActiv, including a complete step-by-step
workflow and a function manual, is available at
<https://goekelab.github.io/proActiv/>.

Additional data on differential promoters in tissues and cancers from
TCGA, ICGC, GTEx, and PCAWG is available at
<https://jglab.org/data-and-software/>.

### Content

  - [Installation](#installation)
  - [Quick Start](#quick-start)
  - [Creating a Promoter Annotation
    object](#creating-a-promoter-annotation-object)
  - [Complete Analysis Workflow: Analyzing Alternative
    Promoters](#complete-analysis-workflow-analyzing-alternative-promoters)
  - [Limitations](#limitations)
  - [Release History](#release-history)
  - [Reference](#reference)
  - [Contributors](#contributors)

### Installation

proActiv can be installed from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("proActiv")
```

### Quick Start

proActiv estimates promoter activity from RNA-Seq data. Promoter
activity is defined as the total amount of transcription initiated at
each promoter. proActiv takes as input either BAM files or junction
files (TopHat2 or STAR), and a promoter annotation object of the
relevant genome. An optional argument `condition` can be supplied,
describing the experimental condition corresponding to each input file.
Here we demonstrate proActiv with STAR junction files (Human genome
GRCh38 GENCODE v34) as input. These files are taken from the [SGNEx
project](https://github.com/GoekeLab/sg-nex-data) but restricted to the
chr1:10,000,000-30,000,000 region, and can be found at
`extdata/vignette`:

  - `extdata/vignette/SGNEx_A549_Illumina_replicate1-run1.subset.SJ.out.tab.gz`
  - `extdata/vignette/SGNEx_A549_Illumina_replicate3-run1.subset.SJ.out.tab.gz`
  - `extdata/vignette/SGNEx_A549_Illumina_replicate5-run1.subset.SJ.out.tab.gz`
  - `extdata/vignette/SGNEx_HepG2_Illumina_replicate2-run1.subset.SJ.out.tab.gz`
  - `extdata/vignette/SGNEx_HepG2_Illumina_replicate4-run1.subset.SJ.out.tab.gz`
  - `extdata/vignette/SGNEx_HepG2_Illumina_replicate5-run1.subset.SJ.out.tab.gz`

<!-- end list -->

``` r
library(proActiv)

## List of STAR junction files as input
files <- list.files(system.file('extdata/vignette/junctions', 
                                package = 'proActiv'), full.names = TRUE)
## Vector describing experimental condition
condition <- rep(c('A549','HepG2'), each=3)
## Promoter annotation for human genome GENCODE v34
promoterAnnotation <- promoterAnnotation.gencode.v34.subset

result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation,
                   condition = condition)
```

`result` is a summarizedExperiment object which can be accessed as
follows:

  - `assays(results)` returns raw/normalized promoter counts,
    absolute/relative promoter activity and gene expression data  
  - `rowData(results)` returns promoter metadata and summarized absolute
    promoter activity by conditions

proActiv can also be run with BAM files as input, but an additional
parameter `genome` must be supplied:

``` r
## From BAM files - genome parameter must be provided
files <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), full.names = TRUE)
result <- proActiv(files = files, 
                   promoterAnnotation = promoterAnnotation.gencode.v34.subset,
                   genome = 'hg38')
```

### Creating a Promoter Annotation object

In order to quantify promoter activity, proActiv uses a set of promoters
based on genome annotations. proActiv allows the creation of a promoter
annotation object for any genome from a TxDb object or from a GTF file
with the `preparePromoterAnnotation` function. Users have the option to
either pass the file path of the GTF/GFF or TxDb to be used, or use the
TxDb object directly as input. proActiv includes pre-calculated promoter
annotations for the human genome (GENCODE v34). However, due to size
constraints, the annotation is restricted to the
chr1:10,000,000-30,000,000 region. Users can build full annotations by
downloading GTF files from [GENCODE](https://www.gencodegenes.org) page
and following the steps below.

Here, we demonstrate creating the subsetted promoter annotation for the
Human genome (GENCODE v34) with both GTF and TxDb:

``` r
## From GTF file path
gtf.file <- system.file('extdata/vignette/annotation/gencode.v34.annotation.subset.gtf.gz', 
                        package = 'proActiv')
promoterAnnotation.gencode.v34.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')
## From TxDb object
txdb.file <- system.file('extdata/vignette/annotation/gencode.v34.annotation.subset.sqlite', 
                         package = 'proActiv')
txdb <- loadDb(txdb.file)
promoterAnnotation.gencode.v34.subset <- preparePromoterAnnotation(txdb = txdb, 
                                                                   species = 'Homo_sapiens')
```

The `PromoterAnnotation` object has 3 slots:

  - `intronRanges` : Intron ranges, giving the corresponding transcripts
    of each intron
  - `promoterIdMapping` : An ID mapping between transcripts, promoter
    IDs and gene IDs  
  - `promoterCoordinates` : Promoter coordinates (TSS) and internal
    promoter state, along with the 3’ coordinate of the first exon

### Complete Analysis Workflow: Analyzing Alternative Promoters

Most human genes have multiple promoters that control the expression of
distinct isoforms. The use of these alternative promoters enables the
regulation of isoform expression pre-transcriptionally. Importantly,
alternative promoters have been found to be important in a wide number
of cell types and diseases. proActiv includes a workflow to identify and
visualize alternative promoter usage between conditions. This workflow
is described in detail
[here](https://goekelab.github.io/proActiv/articles/proActiv.html).

## Release History

**Release 1.1.18**

Date: 7th April 2021

Changes in version 1.1.18:

  - Gene expression data is now stored in the `assays` of the
    summarizedExperiment object returned by `proActiv` to facilitate
    easier filtering of the summarizedExperiment object. The metadata
    slot is now empty.

  - Plotting promoter activity: Implementation of `boxplotPromoters`
    function to plot boxplots of absolute promoter activity, relative
    promoter activity, and gene expression.

  - Identification of alternative promoters: Implementation of
    `getAlternativePromoters`, used to identify promoters that may
    exhibit alternative usage.

**Release 1.0.0**

Release date: 28th October 2020

Released with Bioconductor 3.12

**Release 0.99.0**

Release date: 21st August 2020

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
    `preparePromoterAnnotations` function. The promoter annotation
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
