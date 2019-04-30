context("Load GTF")
library(rucdr)

txdb <- loadGTF("data/test.gtf")

test_that("str_length is number of characters", {
  expect_type(txdb, "S4")
})

file.remove("data/test.gtf.sqlite")
