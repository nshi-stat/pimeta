n1 <- c(30, 45, 60, 25); m1 <- c(-2.1, -1.4, -3.0, -0.8); s1 <- c(2.0, 1.8, 2.5, 1.6)
n2 <- c(30, 44, 58, 26); m2 <- c(-1.0, -1.1, -1.9, -0.9); s2 <- c(2.2, 1.7, 2.4, 1.5)

test_that("mean difference (unpooled variance)", {
  r <- convert_mean(n1, m1, s1, n2, m2, s2)
  expect_s3_class(r, "data.frame")
  expect_equal(r$y, m1 - m2, tolerance = 1e-12)
  expect_equal(r$se, sqrt(s1^2/n1 + s2^2/n2), tolerance = 1e-12)
  # swapping the groups changes only the sign
  s <- convert_mean(n2, m2, s2, n1, m1, s1)
  expect_equal(s$y, -r$y, tolerance = 1e-12)
  expect_equal(s$se, r$se, tolerance = 1e-12)
})

test_that("mean difference (pooled variance)", {
  r <- convert_mean(n1, m1, s1, n2, m2, s2, pooled = TRUE)
  sp2 <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2)
  expect_equal(r$y, m1 - m2, tolerance = 1e-12)
  expect_equal(r$se, sqrt(sp2*(1/n1 + 1/n2)), tolerance = 1e-12)
  # s1 is used (pimeta development versions used s2 twice)
  r2 <- convert_mean(n1, m1, s1*2, n2, m2, s2, pooled = TRUE)
  expect_true(all(r2$se > r$se))
  # equal SDs and equal sizes: pooled = unpooled
  e <- convert_mean(n1, m1, s1, n1, m2, s1, pooled = TRUE)
  u <- convert_mean(n1, m1, s1, n1, m2, s1)
  expect_equal(e$se, u$se, tolerance = 1e-12)
  # agrees with the two-sample t-test
  set.seed(1)
  a <- rnorm(12, 1); b <- rnorm(15)
  t <- stats::t.test(a, b, var.equal = TRUE)
  p <- convert_mean(12, mean(a), sd(a), 15, mean(b), sd(b), pooled = TRUE)
  expect_equal(p$y/p$se, unname(t$statistic), tolerance = 1e-10)
  w <- stats::t.test(a, b, var.equal = FALSE)
  q <- convert_mean(12, mean(a), sd(a), 15, mean(b), sd(b))
  expect_equal(q$y/q$se, unname(w$statistic), tolerance = 1e-10)
})

test_that("negative means are allowed and the input is checked", {
  expect_silent(convert_mean(10, -5, 1, 10, -3, 1))
  expect_error(convert_mean(n1, m1, -s1, n2, m2, s2), "non-negative")
  expect_error(convert_mean(n1[-1], m1, s1, n2, m2, s2), "same length")
  expect_error(convert_mean(c(0, n1[-1]), m1, s1, n2, m2, s2), "positive")
  expect_error(convert_mean(n1, m1, s1, n2, m2, s2, pooled = NA), "TRUE or FALSE")
  expect_error(convert_mean(n1, m1, s1, n2, m2, s2, pooled = "yes"), "TRUE or FALSE")
  expect_error(convert_mean(1, 0, 1, 1, 0, 1, pooled = TRUE), "greater than 2")
})

test_that("clbp and nfld data can be converted", {
  for (nm in c("clbp", "nfld")) {
    e <- new.env()
    utils::data(list = nm, package = "pimeta", envir = e)
    d <- e[[nm]]
    expect_s3_class(d, "data.frame")
    expect_true(all(c("m1", "s1", "n1", "m2", "s2", "n2") %in% names(d)))
    r <- convert_mean(d$n1, d$m1, d$s1, d$n2, d$m2, d$s2)
    expect_equal(nrow(r), nrow(d))
    expect_true(all(is.finite(r$y)) && all(r$se > 0))
  }
  e <- new.env()
  utils::data(list = c("clbp", "nfld"), package = "pimeta", envir = e)
  expect_equal(c(nrow(e$clbp), nrow(e$nfld)), c(23, 25))
})
