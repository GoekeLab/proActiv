context('Promoter Activity and Gene Expression')
library(proActiv)

promoterAnnotation <- promoterAnnotation.gencode.v19
filesTophat <- list.files(system.file('/extdata/tophat2', package = 'proActiv'), full.names = TRUE, pattern = 'sample1|sample2')
fileLabels <- paste0('sample', 1:2)

tophatJunctionCounts <- calculatePromoterReadCounts(promoterAnnotation, filesTophat, fileLabels, 'tophat')
tophatNormCounts <- normalizePromoterReadCounts(tophatJunctionCounts)
tophatAbsProActiv <- getAbsolutePromoterActivity(tophatNormCounts, promoterAnnotation)
tophatGeneExp <- getGeneExpression(tophatAbsProActiv) 
tophatRelProActiv <- getRelativePromoterActivity(tophatAbsProActiv, tophatGeneExp)

test_that('getAbsolutePromoterActivity returns expected output', {
  
  absolutePromoterActivity <- readRDS(system.file('/extdata/testdata/tophat2', 'absolutePromoterActivity.rds', package = 'proActiv'))
  
  expect_type(absolutePromoterActivity, 'list')
  expect_equal(tophatAbsProActiv, absolutePromoterActivity)
  
})


test_that('getGeneExpression returns expected output', {
  
  geneExpression <- readRDS(system.file('/extdata/testdata/tophat2', 'geneExpression.rds', package = 'proActiv'))
  
  expect_type(geneExpression, 'list')
  expect_equal(tophatGeneExp, geneExpression)

})


test_that('getRelativePromoterActivity returns expected output', {

  relativePromoterActivity <- readRDS(system.file('/extdata/testdata/tophat2', 'relativePromoterActivity.rds', package = 'proActiv'))

  expect_type(relativePromoterActivity, 'list')
  expect_identical(tophatRelProActiv, relativePromoterActivity)

})