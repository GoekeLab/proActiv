context('Calculating Promoter Annotation')
library(proActiv)
library(mockery)

gtfPath <- list.files(system.file('extdata/testdata/promoterAnnotation', package = 'proActiv'), full.names = TRUE, pattern = 'gtf')
txdbPath <- list.files(system.file('extdata/testdata/promoterAnnotation', package = 'proActiv'), full.names = TRUE, pattern = 'sqlite')

txdb <- AnnotationDbi::loadDb(txdbPath)

reducedExonRanges <- readRDS(system.file('extdata/testdata/promoterAnnotation', 'reducedExonRanges.rds', package = 'proActiv'))
annotatedIntronRanges <- readRDS(system.file('extdata/testdata/promoterAnnotation', 'annotatedIntronRanges.rds', package = 'proActiv'))
promoterIdMapping <- readRDS(system.file('extdata/testdata/promoterAnnotation', 'promoterIdMapping.rds', package = 'proActiv'))
promoterCoordinates <- readRDS(system.file('extdata/testdata/promoterAnnotation', 'promoterCoordinates.rds', package = 'proActiv'))

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
  
  suppressWarnings(
    promoterAnnotation <- preparePromoterAnnotation(file = gtfPath, species = 'Homo_sapiens')
  )
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
  expect_equal(intronRanges(promoterAnnotation), annotatedIntronRanges[,c('INTRONID', 'TXNAME')])
  expect_equal(promoterIdMapping(promoterAnnotation), promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')])
  expect_equal(promoterCoordinates(promoterAnnotation)[,c('promoterId', 'internalPromoter')], promoterCoordinates[,-2])
  expect_equal(as.numeric(unlist(promoterCoordinates(promoterAnnotation)$intronId)), 
               as.numeric(unlist(reducedExonRanges$intronId)))
  
})

test_that('preparePromoterAnnotation returns expected output with txdb file path',{
  
  promoterAnnotation <- preparePromoterAnnotation(file = txdbPath, species = 'Homo_sapiens')
  
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
  expect_equal(intronRanges(promoterAnnotation), annotatedIntronRanges[,c('INTRONID', 'TXNAME')])
  expect_equal(promoterIdMapping(promoterAnnotation), promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')])
  expect_equal(promoterCoordinates(promoterAnnotation)[,c('promoterId', 'internalPromoter')], promoterCoordinates[,-2])
  expect_equal(as.numeric(unlist(promoterCoordinates(promoterAnnotation)$intronId)), 
               as.numeric(unlist(reducedExonRanges$intronId)))
  
})

test_that('preparePromoterAnnotation returns expected output with txdb object', {
  
  promoterAnnotation <- preparePromoterAnnotation(txdb = txdb, species = 'Homo_sapiens')
  
  expect_s4_class(intronRanges(promoterAnnotation), 'GRanges')
  expect_type(promoterIdMapping(promoterAnnotation), 'list')
  expect_s4_class(promoterCoordinates(promoterAnnotation), 'GRanges')
  
  expect_equal(intronRanges(promoterAnnotation), annotatedIntronRanges[,c('INTRONID', 'TXNAME')])
  expect_equal(promoterIdMapping(promoterAnnotation), promoterIdMapping[,c('transcriptName', 'promoterId', 'geneId')])
  expect_equal(promoterCoordinates(promoterAnnotation)[,c('promoterId', 'internalPromoter')], promoterCoordinates[,-2])
  expect_equal(as.numeric(unlist(promoterCoordinates(promoterAnnotation)$intronId)), 
               as.numeric(unlist(reducedExonRanges$intronId)))
  
})
