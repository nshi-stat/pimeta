# Parametric bootstrap prediction interval (method = "boot")

sbp <- get_data("sbp")
sbp$sigmak <- sbp$se

test_that("bootstrap: 'table' agrees with 'exact'", {
  skip_on_cran()
  ex <- pimeta::pima(sbp$y, sbp$sigmak, B = 10000, seed = 3141592,
                     sampling = "exact")
  tb <- pimeta::pima(sbp$y, sbp$sigmak, B = 10000, seed = 3141592,
                     sampling = "table")
  # the same uniform random numbers are used, so the samples of tau^2
  # agree draw by draw up to the interpolation error
  ok <- is.finite(ex$rnd) & is.finite(tb$rnd)
  expect_lt(max(abs(ex$rnd[ok] - tb$rnd[ok])), 1e-3)
  expect_equal(c(tb$lpi, tb$upi, tb$lci, tb$uci),
               c(ex$lpi, ex$upi, ex$lci, ex$uci), tolerance = 1e-3)
})

test_that("bootstrap: reproducibility", {
  skip_on_cran()
  a <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1)
  b <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1)
  expect_identical(a$rnd, b$rnd)
  expect_identical(c(a$lpi, a$upi), c(b$lpi, b$upi))
  # set.seed() before pima() (seed = NULL) is respected
  set.seed(2)
  c1 <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000)
  set.seed(2)
  c2 <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000)
  expect_identical(c(c1$lpi, c1$upi), c(c2$lpi, c2$upi))
})

test_that("bootstrap: results do not depend on the number of threads", {
  skip_on_cran()
  a <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1, parallel = FALSE)
  b <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1, parallel = 2)
  expect_identical(a$rnd, b$rnd)
})

test_that("bootstrap: user supplied rnd", {
  skip_on_cran()
  a <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1)
  b <- pimeta::pima(sbp$y, sbp$sigmak, rnd = a$rnd, seed = 1)
  expect_identical(a$rnd, b$rnd)
})

test_that("bootstrap: equivariance to the scale and location of y", {
  skip_on_cran()
  a <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1)
  for (f in c(1e-6, 1e4)) {
    b <- pimeta::pima(sbp$y*f, sbp$sigmak*f, B = 5000, seed = 1)
    expect_equal(c(b$lpi, b$upi, b$lci, b$uci)/f, c(a$lpi, a$upi, a$lci, a$uci),
                 tolerance = 1e-6)
    expect_equal(b$rnd/f^2, a$rnd, tolerance = 1e-6)
  }
  d <- pimeta::pima(sbp$y + 1e4, sbp$sigmak, B = 5000, seed = 1)
  expect_equal(c(d$lpi, d$upi, d$lci, d$uci) - 1e4, c(a$lpi, a$upi, a$lci, a$uci),
               tolerance = 1e-6)
})

test_that("bootstrap: no Inf for K = 3 with large heterogeneity", {
  skip_on_cran()
  res <- pimeta::pima(c(-2, 0, 3), c(0.1, 0.1, 0.1), B = 25000, seed = 1)
  expect_equal(sum(!is.finite(res$rnd)), 0)
  expect_true(all(is.finite(c(res$lpi, res$upi, res$lci, res$uci))))
})

test_that("bootstrap: Inf samples of tau^2 (user supplied) are handled as the limit", {
  skip_on_cran()
  a <- pimeta::pima(sbp$y, sbp$sigmak, B = 5000, seed = 1)
  r <- a$rnd
  r[1:10] <- Inf
  expect_warning(b <- pimeta::pima(sbp$y, sbp$sigmak, rnd = r, seed = 1), "tau\\^2 = Inf")
  expect_true(all(is.finite(c(b$lpi, b$upi, b$lci, b$uci))))
  expect_false(any(is.nan(c(b$lpi, b$upi, b$lci, b$uci))))
})

test_that("bootstrap: identical y", {
  skip_on_cran()
  expect_warning(res <- pimeta::pima(rep(0.3, 6), c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3),
                                     B = 5000, seed = 1), "identical")
  expect_true(all(res$rnd == 0))
  expect_equal(c(res$lpi, res$upi, res$lci, res$uci), rep(0.3, 4))
})
