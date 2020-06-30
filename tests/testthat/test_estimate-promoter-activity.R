context('Promoter Activity and Gene Expression')
library(proActiv)

promoterAnnotation <- promoterAnnotation.gencode.v19
filesTophat <- list.files(system.file('/extdata/testdata/tophat2', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
fileLabels <- paste0('sample', 1:2)

tophatJunctionCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat')
tophatNormCounts <- normalizePromoterReadCounts(tophatJunctionCounts)

test_that('getAbsolutePromoterActivity returns expected output', {
  
  absolutePromoterActivity <- readRDS(system.file('/extdata/testdata/tophat2', 'absolutePromoterActivity.rds', package = 'proActiv'))
  
  expect_type(getAbsolutePromoterActivity(tophatNormCounts, promoterAnnotation), 'list')
  expect_equal(getAbsolutePromoterActivity(tophatNormCounts, promoterAnnotation), absolutePromoterActivity)
  
})


test_that('getGeneExpression returns expected output', {
  
  geneExpression <- readRDS(system.file('/extdata/testdata/tophat2', 'geneExpression.rds', package = 'proActiv'))
  tophatAbsProActiv <- getAbsolutePromoterActivity(tophatNormCounts, promoterAnnotation)

  expect_type(getGeneExpression(tophatAbsProActiv), 'list')
  expect_equal(getGeneExpression(tophatAbsProActiv), geneExpression)

})


test_that('getRelativePromoterActivity returns expected output', {

  relativePromoterActivity <- readRDS(system.file('/extdata/testdata/tophat2', 'relativePromoterActivity.rds', package = 'proActiv'))
  tophatAbsProActiv <- getAbsolutePromoterActivity(tophatNormCounts, promoterAnnotation)
  tophatGeneExp <- getGeneExpression(tophatAbsProActiv) 

  expect_type(getRelativePromoterActivity(tophatAbsProActiv, tophatGeneExp), 'list')
  expect_identical(getRelativePromoterActivity(tophatAbsProActiv, tophatGeneExp), relativePromoterActivity)

})
