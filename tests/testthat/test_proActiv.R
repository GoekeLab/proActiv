context('proActiv Wrapper')
library(proActiv)
library(mockery)

promoterAnnotation <- promoterAnnotation.gencode.v19

test_that('proActiv handles non-existent files',{
  
  expect_error(proActiv('', promoterAnnotation))

  })

test_that('proActiv handles multiple file types', {
  
  mockery::stub(proActiv, 'file.exists', function() TRUE)
  files <- c('sample1.bed', 'sample2.junctions')
  expect_error(proActiv(files, promoterAnnotation))

  })

test_that('proActiv handles invalid file types', {
  
  mockery::stub(proActiv, 'file.exists', function() TRUE)
  files <- c('sample1.txt')
  expect_error(proActiv(files, promoterAnnotation))

  })

test_that('proActiv expects genome argument with BAM input', {
  
  mockery::stub(proActiv, 'file.exists', function() TRUE)
  files <- c('sample1.bam')
  expect_error(proActiv(files, promoterAnnotation))
  
})

test_that('proActiv handles invalid condition argument', {
  
  files <- list.files(system.file('extdata/vignette', package = 'proActiv'), 
                      full.names = TRUE, pattern = 'replicate5')
  expect_warning(proActiv(files = files, 
                           promoterAnnotation = promoterAnnotation.gencode.v34,
                           condition = c('A549', 'HepG2', 'MCF7')))
  
})

test_that('proActiv handles compressed input files', {

  files <- list.files(system.file('extdata/vignette', package = 'proActiv'), 
                      full.names = TRUE, pattern = 'replicate5')
  expect_s4_class(proActiv(files = files, 
                           promoterAnnotation = promoterAnnotation.gencode.v34,
                           condition = c('A549', 'HepG2')), 'SummarizedExperiment')

})

test_that('proActiv returns a Summarized Experiment', {

  ### 1) Test Tophat2 BED file input ###
  filesTophat <- list.files(system.file('extdata/testdata/tophat2', package = 'proActiv'),
                            full.names = TRUE, pattern = 'sample1|sample2')
  result <- proActiv(filesTophat, promoterAnnotation)

  promoterCounts <- readRDS(system.file('extdata/testdata/tophat2/promoterCounts.rds', package = 'proActiv'))
  normalizedPromoterCounts <- readRDS(system.file('extdata/testdata/tophat2/normalizedPromoterCounts.rds', package = 'proActiv'))
  absolutePromoterActivity <- readRDS(system.file('extdata/testdata/tophat2/absolutePromoterActivity.rds', package = 'proActiv'))
  geneExpression <- readRDS(system.file('extdata/testdata/tophat2/geneExpression.rds', package = 'proActiv'))
  relativePromoterActivity <- readRDS(system.file('extdata/testdata/tophat2/relativePromoterActivity.rds', package = 'proActiv'))

  expect_s4_class(result, 'SummarizedExperiment')
  expect_identical(promoterCounts, assays(result)$promoterCounts)
  expect_identical(normalizedPromoterCounts, assays(result)$normalizedPromoterCounts)
  expect_identical(absolutePromoterActivity[,c('sample1', 'sample2')], assays(result)$absolutePromoterActivity)
  expect_identical(geneExpression[c('sample1', 'sample2')], metadata(result)$geneExpression)
  expect_identical(relativePromoterActivity[,c('sample1', 'sample2')], assays(result)$relativePromoterActivity)
  expect_identical(absolutePromoterActivity[,c('promoterId', 'geneId')], data.frame(rowData(result)[,c('promoterId', 'geneId')]))

  ## 2) Test STAR junction file input ###
  filesSTAR <- list.files(system.file('extdata/testdata/star', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
  result <- proActiv(filesSTAR, promoterAnnotation)

  promoterCounts <- readRDS(system.file('extdata/testdata/star/promoterCounts.rds', package = 'proActiv'))
  normalizedPromoterCounts <- readRDS(system.file('extdata/testdata/star/normalizedPromoterCounts.rds', package = 'proActiv'))
  absolutePromoterActivity <- readRDS(system.file('extdata/testdata/star/absolutePromoterActivity.rds', package = 'proActiv'))
  geneExpression <- readRDS(system.file('extdata/testdata/star/geneExpression.rds', package = 'proActiv'))
  relativePromoterActivity <- readRDS(system.file('extdata/testdata/star/relativePromoterActivity.rds', package = 'proActiv'))

  expect_s4_class(result, 'SummarizedExperiment')
  expect_identical(promoterCounts, assays(result)$promoterCounts)
  expect_identical(normalizedPromoterCounts, assays(result)$normalizedPromoterCounts)
  expect_identical(absolutePromoterActivity[,c('sample1', 'sample2')], assays(result)$absolutePromoterActivity)
  expect_identical(geneExpression[c('sample1', 'sample2')], metadata(result)$geneExpression)
  expect_identical(relativePromoterActivity[,c('sample1', 'sample2')], assays(result)$relativePromoterActivity)
  expect_identical(absolutePromoterActivity[,c('promoterId', 'geneId')], data.frame(rowData(result)[,c('promoterId', 'geneId')]))

  # 3) Test BAM file input
  bamfiles <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), full.names = TRUE)
  suppressWarnings(
   result <- proActiv(bamfiles, promoterAnnotation.gencode.v34, genome = 'hg38')
  )
  expect_s4_class(result, 'SummarizedExperiment')
  expect_identical(dim(result)[1], length(promoterCoordinates(promoterAnnotation.gencode.v34)))

})
