context('Promoter and Junction Read Counts')
library(proActiv)
library(mockery)

promoterAnnotation <- promoterAnnotation.gencode.v19
promoterCoordinates <- promoterCoordinates(promoterAnnotation)
intronRanges <- intronRanges(promoterAnnotation)

filesTophat <- list.files(system.file('/extdata/tophat2', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
filesSTAR <- list.files(system.file('/extdata/star', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
fileLabels <- paste0('sample', 1:2)

test_that('calculateJunctionReadCounts returns expected output',{

  ### Tophat input ###
  promoterCounts <- readRDS(system.file('/extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))
  tophatJunctionCounts <- calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesTophat[1], 'tophat')

  expect_type(tophatJunctionCounts, 'double')
  expect_equal(as.numeric(tophatJunctionCounts), promoterCounts$sample1)

  ### STAR input ###
  promoterCounts <- readRDS(system.file('/extdata/testdata/star', 'promoterCounts.rds', package = 'proActiv'))
  starJunctionCounts <- calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesSTAR[1], 'star')

  expect_type(starJunctionCounts, 'double')
  expect_equal(as.numeric(starJunctionCounts), promoterCounts$sample1)

  ### BAM input ###

})

tophatPromoterCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat')
starPromoterCounts <- calculatePromoterReadCounts(promoterAnnotation, filesSTAR , fileLabels, 'star')

test_that('calculatePromoterReadCounts returns expected output', {

  ### 1) Tophat input ###
  promoterCounts <- readRDS(system.file('/extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(tophatPromoterCounts, 'list')
  expect_equal(tophatPromoterCounts, promoterCounts)

  ### 2) STAR input ###
  promoterCounts <- readRDS(system.file('/extdata/testdata/star', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(starPromoterCounts, 'list')
  expect_equal(starPromoterCounts, promoterCounts)

  ### 3) BAM input ###

})

test_that('parallelised calculatePromoterReadCounts returns correct output', {

  promoterCounts <- readRDS(system.file('/extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))

  suppressWarnings(
    tophatPromoterCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat', numberOfCores = 2)
  )
  expect_equal(tophatPromoterCounts, promoterCounts)

})

test_that('normalizedReadCounts returns expected output', {

  normalizedPromoterCounts <- readRDS(system.file('/extdata/testdata/tophat2', 'normalizedPromoterCounts.rds', package = 'proActiv'))
  normalizedCounts <- normalizePromoterReadCounts(tophatPromoterCounts)

  expect_type(normalizedCounts, 'list')
  expect_equal(normalizedCounts, normalizedPromoterCounts)

})

test_that('normalizedReadCounts throws error when DESeq2 unavailable', {

  mockery::stub(normalizePromoterReadCounts, 'requireNamespace', function() FALSE)
  expect_error(normalizePromoterReadCounts(tophatPromoterCounts))

})