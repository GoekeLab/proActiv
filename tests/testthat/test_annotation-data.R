context('Calculating Promoter Annotation')
library(proActiv)
library(mockery)

gtfPath <- system.file('extdata/vignette/annotations/gencode.v34.annotation.subset.gtf.gz', package = 'proActiv')
txdbPath <- system.file('extdata/vignette/annotations/gencode.v34.annotation.subset.sqlite', package = 'proActiv')

txdb <- loadDb(txdbPath)

test_that('preparePromoterAnnotation handles ambiguous argument specification', {
  expect_error(preparePromoterAnnotation(gtfPath, 'Homo_sapiens'))
})

test_that('preparePromoterAnnotation handles multiple file types', {
  
  mockery::stub(preparePromoterAnnotation, 'file.exists', function() TRUE)
  files <- 'sample.gtf'
  expect_error(preparePromoterAnnotation(file = file, species = 'Homo_sapiens'))
  
})

test_that('preparePromoterAnnotation handles non-existent files', {
  
  file <- 'sample.txt'
  expect_error(preparePromoterAnnotation(file = file, species = 'Homo_sapiens'))
  
})

test_that('preparePromoterAnnotation handles invalid file types', {
  
  mockery::stub(preparePromoterAnnotation, 'file.exists', function() TRUE)
  file <- 'sample.txt'
  expect_error(preparePromoterAnnotation(file = file, species = 'Homo_sapiens'))
  
})


test_that('preparePromoterAnnotation returns expected output with gtf file path',{
  
  ## Windows OS throws warnings
  suppressWarnings(
    promoterAnnotation <- preparePromoterAnnotation(file = gtfPath, species = 'Homo_sapiens')
  )
  expect_identical(promoterAnnotation.gencode.v34.subset, promoterAnnotation)
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
})

test_that('preparePromoterAnnotation returns expected output with txdb file path',{
  
  promoterAnnotation <- preparePromoterAnnotation(file = txdbPath, species = 'Homo_sapiens')
  expect_identical(promoterAnnotation.gencode.v34.subset, promoterAnnotation)
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
})

test_that('preparePromoterAnnotation returns expected output with txdb object', {
  
  promoterAnnotation <- preparePromoterAnnotation(txdb = txdb, species = 'Homo_sapiens')
  expect_identical(promoterAnnotation.gencode.v34.subset, promoterAnnotation)
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
})
