context("filter")
library(rucdr)
library(magrittr)
library(dplyr)

pl <- pipeline()

pl$metadata <- data.frame(id=c("A","B","C","D"))

pl  <- pl %>% filter(id != "B")

test_that("salmon()", {
    expect_s3_class(pl$metadata.selection, "data.frame")
    expect_true(nrow(pl$metadata.selection) == 3)
    expect_true(! "B" %in%  pl$metadata.selection$id)
})
