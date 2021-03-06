% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proActiv.R
\name{proActiv}
\alias{proActiv}
\title{Estimates promoter counts and activity in a single command}
\usage{
proActiv(
  files,
  promoterAnnotation,
  fileLabels = NULL,
  condition = NULL,
  genome = NULL,
  ncores = 1
)
}
\arguments{
\item{files}{A character vector. The list of input files for 
which the junction read counts will be calculated}

\item{promoterAnnotation}{A PromoterAnnotation object containing the
intron ranges, promoter coordinates and the promoter id mapping}

\item{fileLabels}{A character vector. The labels of input files 
for which the junction read counts will be calculated. These labels will be 
used as column names for each output data.frame object. If not provided,
filenames will be used as labels. Defaults to NULL}

\item{condition}{A character vector. The condition to which each sample
belong to. Must correspond to the order of the files. If supplied, 
results are summarized by condition. Defaults to NULL}

\item{genome}{A character. Genome version. Must be specified if input file
type is a BAM file. Defaults to NULL}

\item{ncores}{A numeric value. The number of cores to be used for 
counting junction reads. Defaults to 1 (no parallelization). This parameter 
will be used as an argument to BiocParallel::bplapply}
}
\value{
A SummarizedExperiment object with assays giving promoter counts 
  and activity with gene expression. rowData contains
  promoter metadata and absolute promoter activity summarized across
  conditions (if condition is provided)
}
\description{
Estimates promoter counts and activity in a single command
}
\examples{

files <- list.files(system.file('extdata/vignette/junctions', 
                       package = 'proActiv'), 
                       full.names = TRUE, pattern = 'replicate5')
promoterAnnotation <- promoterAnnotation.gencode.v34.subset
result <- proActiv(files = files,
                       promoterAnnotation  = promoterAnnotation,
                       condition = rep(c('A549', 'HepG2'), each=1),
                       fileLabels = NULL,
                       ncores = 1)
                           
}
