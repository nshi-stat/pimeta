d <- get_data("sbp")

test_that("invalid inputs are rejected", {
  expect_error(pima(d$y), "Either 'se' or 'v'")
  expect_error(pima(c(d$y[-1], NA), d$se, method = "HTS"), "missing")
  expect_error(pima(c(d$y[-1], NaN), d$se, method = "HTS"), "NaN")
  expect_error(pima(c(d$y[-1], Inf), d$se, method = "HTS"), "infinite")
  expect_error(pima(d$y, -d$se, method = "HTS"), "non-negative")
  expect_error(pima(d$y, d$se[-1], method = "HTS"), "same length")
  expect_error(pima(d$y, d$se, alpha = 1.5, method = "HTS"), "range")
  expect_error(pima(d$y, d$se, method = "XX"))
  expect_error(pima(d$y, d$se, lower = 2, upper = 1), "greater than 'lower'")
  expect_error(cima(d$y, d$se[-1], method = "DL"), "same length")
  expect_error(cima(d$y, d$se, method = "XX"))
  expect_error(tau2h(d$y, d$se[-1]), "same length")
  expect_error(tau2h(d$y, d$se, method = "XX"))
  expect_error(convert_bin(1:2, c(10, 10), 1:2, 10), "same length")
  expect_error(convert_bin(-1, 10, 1, 10), "non-negative")
})

test_that("'v' is equivalent to 'se'", {
  for (m in c("HTS", "HK", "SJ", "KR", "APX")) {
    expect_equal(pima(d$y, v = d$se^2, method = m)[c("lpi", "upi", "lci", "uci")],
                 pima(d$y, d$se, method = m)[c("lpi", "upi", "lci", "uci")],
                 tolerance = 1e-12)
  }
  expect_equal(cima(d$y, v = d$se^2, method = "PL")[c("lci", "uci")],
               cima(d$y, d$se, method = "PL")[c("lci", "uci")],
               tolerance = 1e-12)
})
