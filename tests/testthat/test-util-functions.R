library(dispersim)
context("Utility functions")

test_that("isPositiveInteger does its job", {
  expect_equal(isPositiveInteger(1), TRUE)
  expect_equal(isPositiveInteger(1000), TRUE)
  expect_equal(isPositiveInteger(0), FALSE)
  expect_equal(isPositiveInteger(0.000001), FALSE)
  expect_equal(isPositiveInteger(-1), FALSE)
})


test_that("checkParams catches errors", {
  

})

