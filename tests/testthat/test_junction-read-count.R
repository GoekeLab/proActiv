context('Promoter and Junction Read Counts')
library(proActiv)
library(mockery)

promoterAnnotation <- promoterAnnotation.gencode.v19
promoterCoordinates <- promoterCoordinates(promoterAnnotation)
intronRanges <- intronRanges(promoterAnnotation)

filesTophat <- list.files(system.file('extdata/testdata/tophat2', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
filesSTAR <- list.files(system.file('extdata/testdata/star', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
bamfiles <- list.files(system.file('extdata/testdata/bam', package = 'proActiv'), full.names = TRUE)
fileLabels <- paste0('sample', 1:2)

test_that('calculateJunctionReadCounts returns expected output',{

  ### Tophat input ###
  promoterCounts <- readRDS(system.file('extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesTophat[1], 'tophat'), 'double')
  expect_equal(as.numeric(calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesTophat[1], 'tophat')), promoterCounts$sample1)

  ### STAR input ###
  promoterCounts <- readRDS(system.file('extdata/testdata/star', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesSTAR[1], 'star'), 'double')
  expect_equal(as.numeric(calculateJunctionReadCounts(promoterCoordinates, intronRanges, filesSTAR[1], 'star')), promoterCounts$sample1)

  ### BAM input ###
  expect_type(
    suppressWarnings(calculateJunctionReadCounts(promoterCoordinates(promoterAnnotation.gencode.v34), 
                                          intronRanges(promoterAnnotation.gencode.v34),
                                          bamfiles[1], 'bam', 'hg38')), 'double')
  expect_equal(length(
    suppressWarnings(calculateJunctionReadCounts(promoterCoordinates(promoterAnnotation.gencode.v34),
                                                 intronRanges(promoterAnnotation.gencode.v34),
                                                 bamfiles[1], 'bam', 'hg38'))), 
    length(promoterCoordinates(promoterAnnotation.gencode.v34)))

})


test_that('calculatePromoterReadCounts returns expected output', {

  ### 1) Tophat input ###
  promoterCounts <- readRDS(system.file('extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat'), 'list')
  expect_equal(calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat'), promoterCounts)

  ### 2) STAR input ###
  promoterCounts <- readRDS(system.file('extdata/testdata/star', 'promoterCounts.rds', package = 'proActiv'))

  expect_type(calculatePromoterReadCounts(promoterAnnotation, filesSTAR , fileLabels, 'star'), 'list')
  expect_equal(calculatePromoterReadCounts(promoterAnnotation, filesSTAR , fileLabels, 'star'), promoterCounts)

  ### 3) BAM input ###
  expect_type(
    suppressWarnings(calculatePromoterReadCounts(promoterAnnotation.gencode.v34, bamfiles, fileLabels, 'bam', 'hg38')), 'list'
  )
  expect_equal(
    suppressWarnings(dim(
      calculatePromoterReadCounts(promoterAnnotation.gencode.v34, bamfiles, fileLabels, 'bam', 'hg38'))), 
    c(length(promoterCoordinates(promoterAnnotation.gencode.v34)),2))

})

test_that('parallelised calculatePromoterReadCounts returns correct output', {

  promoterCounts <- readRDS(system.file('extdata/testdata/tophat2', 'promoterCounts.rds', package = 'proActiv'))

  suppressWarnings(
    expect_equal(calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat', numberOfCores = 2), promoterCounts)
  )

})

test_that('normalizedReadCounts returns expected output', {

  normalizedPromoterCounts <- readRDS(system.file('extdata/testdata/tophat2', 'normalizedPromoterCounts.rds', package = 'proActiv'))
  tophatPromoterCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat')

  expect_type(normalizePromoterReadCounts(tophatPromoterCounts), 'list')
  expect_equal(normalizePromoterReadCounts(tophatPromoterCounts), normalizedPromoterCounts)

})

test_that('normalizedReadCounts throws error when DESeq2 unavailable', {

  mockery::stub(normalizePromoterReadCounts, 'requireNamespace', function() FALSE)
  tophatPromoterCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat')
  expect_error(normalizePromoterReadCounts(tophatPromoterCounts))

})
