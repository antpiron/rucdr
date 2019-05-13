context("Pipeline")
library(rucdr)
library(magrittr)


pl <- data.frame(id=character()) %>% pipeline()

test_that("pl is a pipeline()", {
    expect_s3_class(pl, "pipeline")
    expect_true(! is.null(pl$option))
})

