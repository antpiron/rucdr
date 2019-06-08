context("filter")
library(rucdr)
library(magrittr)


test_that("filter()", {
    
    metadata <- data.frame(id=c("A","B","C","D"),
                              quant.sf.fn=c("A","B","C","D"))
    
    row.names(metadata) <- metadata$id

    pl <- metadata %>% pipeline()
    
    pl  <- pl %>% filter(id != "B")
    
    expect_s3_class(pl$metadata.filtered, "data.frame")
    expect_true(nrow(pl$metadata.filtered) == 3)
    expect_true(! "B" %in%  pl$metadata.filtered$id)
    expect_true(all(pl$metadata.filtered$id ==
                    row.names(pl$metadata.filtered)))

    pl  <- pl %>% filter(id != "A", unfiltered=T)

    expect_true("B" %in%  pl$metadata.filtered$id)
    expect_true(! "A" %in%  pl$metadata.filtered$id)
})
