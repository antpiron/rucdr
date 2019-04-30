context("Expression heatmap")
library(rucdr)

test_that("expression.heatmap is a list", {
  expect_type(expression.heatmap(matrix(rnorm(64), nrow=8)), "list")
})

