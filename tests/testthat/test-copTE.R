test_that("copTE works", {

  Z <- matrix(rnorm(2000),ncol = 2)
  res <- copTE(x = Z[,1],y = Z[,2],lx = 1,ly = 1,nboot = 0)

  expect_equal(length(res), 2)
  expect_equal(length(res$estimate), 2)
  expect_equal(length(res$pval), 2)
})
