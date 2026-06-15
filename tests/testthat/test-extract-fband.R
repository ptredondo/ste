test_that("extract.fband works", {
  series <- rnorm(1000)
  series.fband <- extract.fband(x = series,samprate = 100,fband = c("delta","theta"))
  expect_equal(length(series.fband), 2)
  expect_equal(length(series.fband[[1]]), 1000)
  expect_equal(length(series.fband[[2]]), 1000)
})
