m1 <- c(10, 0, 25); n1 <- c(50, 30, 25)
m2 <- c(5, 3, 20); n2 <- c(50, 30, 40)

test_that("log odds ratio (with 0.5 continuity correction)", {
  r <- convert_bin(m1, n1, m2, n2, type = "logOR")
  a <- m1 + 0.5; b <- n1 - m1 + 0.5; c <- m2 + 0.5; d <- n2 - m2 + 0.5
  expect_equal(r$y, log(a*d/(b*c)), tolerance = 1e-12)
  expect_equal(r$se, sqrt(1/a + 1/b + 1/c + 1/d), tolerance = 1e-12)
  # swapping the groups changes only the sign
  s <- convert_bin(m2, n2, m1, n1, type = "logOR")
  expect_equal(s$y, -r$y, tolerance = 1e-12)
  expect_equal(s$se, r$se, tolerance = 1e-12)
})

test_that("log relative risk", {
  r <- convert_bin(m1, n1, m2, n2, type = "logRR")
  expect_equal(r$y, log((m1 + 0.5)/(n1 + 0.5)) - log((m2 + 0.5)/(n2 + 0.5)),
               tolerance = 1e-12)
  expect_true(all(r$se > 0))
  s <- convert_bin(m2, n2, m1, n1, type = "logRR")
  expect_equal(s$y, -r$y, tolerance = 1e-12)
})

test_that("risk difference", {
  r <- convert_bin(m1, n1, m2, n2, type = "RD")
  expect_equal(r$y, m1/n1 - m2/n2, tolerance = 1e-12)
  expect_true(all(r$se > 0))
  expect_s3_class(r, "data.frame")
})
