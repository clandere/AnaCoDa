library(testthat)
library(AnaCoDa)

context("Parameter")

test_that("general parameter functions", {
  expect_equal(testParameter(), 0)
})

# TODO: Implement the following
p <-new(Parameter)
