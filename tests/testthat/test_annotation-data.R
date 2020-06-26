context('Calculating Promoter Annotation')
library(proActiv)
library(mockery)

gtf <- list.files(system.file('/extdata/testdata/promoterAnnotation', package = 'proActiv'), full.names = TRUE, pattern = 'gtf')
txdb <- list.files(system.file('/extdata/testdata/promoterAnnotation', package = 'proActiv'), full.names = TRUE, pattern = 'sqlite')

reducedExonRanges <- readRDS(system.file('/extdata/testdata/promoterAnnotation', 'reducedExonRanges.rds', package = 'proActiv'))
annotatedIntronRanges <- readRDS(system.file('/extdata/testdata/promoterAnnotation', 'annotatedIntronRanges.rds', package = 'proActiv'))
promoterIdMapping <- readRDS(system.file('/extdata/testdata/promoterAnnotation', 'promoterIdMapping.rds', package = 'proActiv'))
promoterCoordinates <- readRDS(system.file('/extdata/testdata/promoterAnnotation', 'promoterCoordinates.rds', package = 'proActiv'))

test_that('proActiv handles multiple file types', {
  
  mockery::stub(proActiv, 'file.exists', function() TRUE)
  files <- c('sample1.bed', 'sample2.junctions')
  expect_error(proActiv(files, promoterAnnotation))
  
})

test_that('preparePromoterAnnotation handles non-existent files', {
  
  file <- 'tmp.txt'
  expect_error(preparePromoterAnnotation(file, 'Homo_sapiens'))
  
})

test_that('preparePromoterAnnotation handles invalid file types', {
  
  mockery::stub(proActiv, 'file.exists', function() TRUE)
  file <- c('tmp.txt')
  expect_error(preparePromoterAnnotation(file, 'Homo_sapiens'))
  
})


test_that('preparePromoterAnnotation returns expected output with gtf',{
  
  suppressWarnings(
    promoterAnnotation <- preparePromoterAnnotation(gtf, 'Homo_sapiens')
  )
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
  expect_equal(intronRanges(promoterAnnotation), annotatedIntronRanges[,c('INTRONID', 'TXNAME')])
  expect_equal(promoterIdMapping(promoterAnnotation), promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')])
  expect_equal(promoterCoordinates(promoterAnnotation)[,c('promoterId', 'geneId', 'internalPromoter')], promoterCoordinates)
  expect_equal(mcols(promoterCoordinates(promoterAnnotation))$intronId, mcols(reducedExonRanges)$intronId)
  
})

test_that('preparePromoterAnnotation returns expected output with txdb',{
  
  suppressWarnings(
    promoterAnnotation <- preparePromoterAnnotation(txdb, 'Homo_sapiens')
  )
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
  expect_equal(intronRanges(promoterAnnotation), annotatedIntronRanges[,c('INTRONID', 'TXNAME')])
  expect_equal(promoterIdMapping(promoterAnnotation), promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')])
  expect_equal(promoterCoordinates(promoterAnnotation)[,c('promoterId', 'geneId', 'internalPromoter')], promoterCoordinates)
  expect_equal(mcols(promoterCoordinates(promoterAnnotation))$intronId, mcols(reducedExonRanges)$intronId)
  
})
