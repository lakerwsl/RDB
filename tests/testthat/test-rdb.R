context("Testing rdb function")
library(RDB)
library(testthat)




test_that("'rdb' function provides expected results", {
  set.seed(1)
  m=50
  d=100 
  P=matrix(runif(m*d),nrow=m,ncol=d)
  Z=rep(0,m)
  Z[1:(m/2)]=1
  res = rdb(P,Z)
  expect_equal(is.null(res), F)
})
