context('proActiv Wrapper')
library(proActiv)

promoterAnnotation <- promoterAnnotation.gencode.v34.subset

test_that('proActiv handles non-existent files',{
  
  expect_error(proActiv('', promoterAnnotation))

  })

test_that('proActiv expects genome argument with BAM input', {
  
  bam <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'),
                         full.names = TRUE)
  expect_error(proActiv(bam, promoterAnnotation))
  
})

test_that('proActiv handles invalid condition argument', {
  
  files <- list.files(system.file('extdata/vignette/junctions', package = 'proActiv'), 
                      full.names = TRUE, pattern = 'replicate5')
  expect_warning(proActiv(files = files, 
                           promoterAnnotation = promoterAnnotation,
                           condition = c('A549', 'HepG2', 'MCF7')))
  
})

test_that('proActiv returns a Summarized Experiment with junction input', {

  ## Test junction file input ###
  junctions <- list.files(system.file('extdata/vignette/junctions', package = 'proActiv'), full.names = TRUE)
  expect_s4_class(proActiv(junctions, promoterAnnotation), 'SummarizedExperiment')
  expect_s4_class(proActiv(junctions, promoterAnnotation, condition = rep(c('A549', 'HepG2'), each=3)), 'SummarizedExperiment')
  
  result <- proActiv(junctions, promoterAnnotation, condition = rep(c('A549', 'HepG2'), each=3))  
  expect_equal(length(assays(result)), 4)
  expect_equal(ncol(rowData(result)), 12)
  
})

# test_that('proActiv returns a Summarized Experiment with bam input', {
#   
#   ## Test BAM file input
#   bams <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), full.names = TRUE)
#   suppressWarnings(
#       expect_s4_class(proActiv(bams, promoterAnnotation, genome = 'hg38'), 'SummarizedExperiment')
#   )
#   suppressWarnings(
#       result <- proActiv(bams, promoterAnnotation, genome = 'hg38')
#   )
#   expect_equal(length(assays(result)), 4)
#   expect_equal(ncol(rowData(result)), 8)
#   
# })


test_that('proActiv parallelisation returns expected output', {
  
  junctions <- list.files(system.file('extdata/vignette/junctions', package = 'proActiv'), full.names = TRUE)
  ## Windows OS throws warnings
  suppressWarnings(
      expect_identical(proActiv(junctions, promoterAnnotation, ncores = 2), proActiv(junctions, promoterAnnotation, ncores = 1))
  )
  
})
