test_that("causal.fband works", {

  phi_tg <- matrix(c(-0.1,    0,    0,    0,
                     2.5,  0.1,    0,    0,
                     0,    0, -0.1,    0,
                     2.5,    0,  2.5,  0.1), 4, byrow = TRUE)
  Z <- causal.fband(nsec = 10,samprate = 100,fband = c("theta","theta","gamma","gamma"),
                    phi = phi_tg,sigma = diag(rep(1,4)),range = c(20,30))
  expect_equal(nrow(Z), 1000)
  expect_equal(ncol(Z), 4)
})
