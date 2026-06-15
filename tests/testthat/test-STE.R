test_that("STE works", {

  Z <- matrix(rnorm(5000),ncol = 2)
  res <- STE(x = Z[,1],y = Z[,2],samprate = 100,blength = 50,
           fbandx = "theta",fbandy = "theta",lx = 1,ly = 1,nboot = 0)
  expect_equal(length(res), 2)
  expect_equal(length(res$estimate), 2)
  expect_equal(length(res$pval), 2)

})
