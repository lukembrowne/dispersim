library(dispersim)
context("Utility functions")

test_that(".isPositiveInteger does its job", {
  expect_equal(.isPositiveInteger(1), TRUE)
  expect_equal(.isPositiveInteger(1000), TRUE)
  expect_equal(.isPositiveInteger(0), FALSE)
  expect_equal(.isPositiveInteger(0.000001), FALSE)
  expect_equal(.isPositiveInteger(-1), FALSE)
})


test_that(".checkParams catches errors", {
  

})

test_that(".pasteLoci does its job..", {
  expect_equal(.pasteLoci(cbind(1,1)), matrix("1-1"))
  expect_equal(.pasteLoci(cbind(c(1,2),c(1,2))), matrix(c("1-1", "2-2")))
  expect_equal(.pasteLoci(cbind(c(1,2),c(1,2), c(3,4), c(3,4))), 
               matrix(c("1-1", "2-2", "3-3", "4-4"), nrow = 2))
})
