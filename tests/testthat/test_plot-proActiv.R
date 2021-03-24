context('plotPromoters')
library(proActiv)

promoterAnnotation <- promoterAnnotation.gencode.v34.subset

files <- list.files(system.file('extdata/vignette/junctions',
                                package = 'proActiv'), full.names = TRUE)

result <- proActiv(files = files,
                    promoterAnnotation  = promoterAnnotation,
                    condition = rep(c('A549','HepG2'), each=3))

txdb <- loadDb(system.file('extdata/vignette/annotations',
                            'gencode.v34.annotation.rap1gap.sqlite',
                            package = 'proActiv'))

ranges <- readRDS(system.file('extdata/vignette/annotations',
                              'exonsBy.rap1gap.rds',
                              package = 'proActiv'))
gene <- 'ENSG00000076864.19'

test_that('plotPromoters throws error if gene not found', {
    
    expect_error(plotPromoters(result, gene = '!!??', txdb = txdb))
    
})

test_that('plotPromoters throws error if both txdb and ranges provided', {
    
    expect_error(plotPromoters(result, gene, txdb = txdb, ranges = ranges))
    
}) 
