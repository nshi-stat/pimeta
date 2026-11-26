test_that("pwchisq reduces to pchisq for a single central chi-square", {
  for (df in 1:4) {
    for (x in c(0.5, 2, 5, 12)) {
      expect_equal(pwchisq(x, lambda = 1, nu = df), stats::pchisq(x, df),
                   tolerance = 1e-6)
      expect_equal(pwchisq(x, lambda = 2.5, nu = df), stats::pchisq(x/2.5, df),
                   tolerance = 1e-6)
    }
  }
})

test_that("pwchisq for a sum of chi-squares with equal weights", {
  for (x in c(1, 4, 9)) {
    expect_equal(pwchisq(x, lambda = c(1, 1), nu = c(1, 2), delta = c(0, 0)),
                 stats::pchisq(x, 3), tolerance = 1e-6)
  }
})

test_that("pwchisq for a linear combination agrees with numerical convolution", {
  # P(X1 + 2 X2 <= x), X1, X2 ~ chi^2_2
  conv <- function(x) {
    stats::integrate(function(u) stats::dchisq(u, 2)*stats::pchisq((x - u)/2, 2),
                     0, x, rel.tol = 1e-10)$value
  }
  for (x in c(1, 5, 15)) {
    expect_equal(pwchisq(x, lambda = c(1, 2), nu = c(2, 2), delta = c(0, 0)),
                 conv(x), tolerance = 1e-6)
  }
})

test_that("pwchisq reproduces Table 1 (Q6) of Farebrother (1984)", {
  p <- sapply(c(20, 100, 200), function(x)
    pwchisq(x, lambda = c(7, 3), nu = c(6, 2), delta = c(6, 2)))
  expect_equal(round(p, 4), c(0.0061, 0.5913, 0.9779))
})

test_that("pwchisq checks the lengths", {
  expect_error(pwchisq(1, lambda = c(1, 2), nu = 1, delta = 0), "same length")
})
