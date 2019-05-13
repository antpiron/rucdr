context("Pipeline options")
library(rucdr)
library(magrittr)


pl <- data.frame(id=character()) %>%
    pipeline() %>% options(nthreads=1000,
                           plop="Plop!")

test_that("pl is a pipeline()", {
    expect_s3_class(pl, "pipeline")
    expect_true(! is.null(pl$option))
    expect_equal(pl$option$nthreads, 1000)
    expect_equal(pl$option$plop, "Plop!")
})

