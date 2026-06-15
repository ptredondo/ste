test_that("sim.fband works", {

  series.fband <- sim.fband(nsec = 10,samprate = 100,fband = c("delta","theta"))
  expect_equal(nrow(series.fband), 1000)
  expect_equal(ncol(series.fband), 2)
})
