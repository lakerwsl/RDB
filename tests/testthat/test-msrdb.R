context("Testing MsRDB function")
library(RDB)
library(testthat)

data(vangay)


test_that("'MsRDB' function provides expected results", {
  res = msrdb(seqtabInterest,Z, knn = 20)
  expect_equal(is.null(res), F)
})
