context("Expression heatmap")
library(rucdr)

test_that("str_length is number of characters", {
  expect_type(expression.heatmap(matrix(rnorm(64), nrow=8)), "list")
})

