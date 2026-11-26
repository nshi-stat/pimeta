test_that("I^2 with the DL estimator equals (Q - (K - 1))/Q", {
  for (d in datasets()) {
    k <- length(d$y)
    q <- ref_q(d$y, d$se^2)
    t <- tau2h(d$y, d$se, method = "DL")$tau2h
    expected <- if (q > k - 1) 100*(q - (k - 1))/q else 0
    expect_equal(i2h(d$se, t), expected, tolerance = 1e-10)
  }
})

test_that("I^2 agrees with metafor", {
  skip_if_not_installed("metafor")
  for (d in datasets()) {
    fit <- metafor::rma(yi = d$y, sei = d$se, method = "REML")
    expect_equal(i2h(d$se, fit$tau2), fit$I2, tolerance = 1e-8)
  }
})

test_that("I^2 is in [0, 100) and increasing in tau^2", {
  se <- get_data("sbp")$se
  i <- sapply(c(0, 0.01, 0.1, 1, 10), function(t) i2h(se, t))
  expect_equal(i[1], 0)
  expect_true(all(diff(i) > 0))
  expect_lt(max(i), 100)
})
