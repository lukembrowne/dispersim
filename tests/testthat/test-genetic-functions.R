library(dispersim)
context("Genetic calculations")

test_that("calcAlleleFreq returns proper values", {
  expect_equal(as.numeric(calcAlleleFreq(c(1,1,1))), 1)
  expect_equal(as.numeric(calcAlleleFreq(c(1,1,2,2))), c(0.5, 0.5))
  expect_equal(as.numeric(calcAlleleFreq(c(1,1,2,2)), c(1,1,2,2)), c(0.5, 0.5))
  expect_equal(as.numeric(calcAlleleFreq(c(1,1,2,2, NA), c(1,1,2,2, NA))), c(0.5, 0.5))
})


test_that("calcHe returns proper values", {
  expect_equal(calcHe(c(1,1,1)), 0)
  expect_equal(calcHe(c(1,1,1, NA)), 0)
  expect_equal(calcHe(c(1,2,3)), 1)
  expect_equal(calcHe(c(1,1,2)), (2/3))
  expect_equal(calcHe(c(1,1), c(1,1)), 0)
  expect_equal(calcHe(c(1,2), c(3,4)), 1)
  expect_equal(calcHe(c(1,1), c(2,2)), (2/3))
})