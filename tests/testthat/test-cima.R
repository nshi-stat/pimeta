test_that("Wald-type CIs of cima agree with those of pima", {
  for (d in datasets()) {
    pairs <- c(DL = "HTS", HK = "HK", SJ = "SJ", KR = "KR", APX = "APX")
    for (m in names(pairs)) {
      c1 <- cima(d$y, d$se, method = m)
      p1 <- pima(d$y, d$se, method = pairs[[m]])
      expect_equal(c(c1$muhat, c1$lci, c1$uci), c(p1$muhat, p1$lci, p1$uci),
                   tolerance = 1e-12)
      expect_equal(c1$K, length(d$y))
    }
  }
})

test_that("profile likelihood CI satisfies the likelihood ratio condition", {
  crit <- stats::qchisq(0.95, 1)/2
  for (d in datasets()) {
    v <- d$se^2
    r <- cima(d$y, d$se, method = "PL")
    expect_close(r$tau2h, ref_tau2(d$y, d$se, "ML"), tol = 1e-4)
    lmax <- ref_pll(r$muhat, d$y, v)
    expect_equal(lmax - ref_pll(r$lci, d$y, v), crit, tolerance = 1e-2)
    expect_equal(lmax - ref_pll(r$uci, d$y, v), crit, tolerance = 1e-2)
    expect_true(r$lci < r$muhat && r$muhat < r$uci)
  }
})

test_that("Bartlett-type corrected CI is wider than the profile likelihood CI", {
  for (d in datasets()) {
    pl <- cima(d$y, d$se, method = "PL")
    bc <- cima(d$y, d$se, method = "BC")
    expect_equal(bc$muhat, pl$muhat)
    expect_lt(bc$lci, pl$lci)
    expect_gt(bc$uci, pl$uci)
  }
})

test_that("bootstrap CI of cima is that of pima (tol and sampling are passed)", {
  skip_on_cran()
  d <- get_data("sbp")
  for (s in c("table", "exact")) {
    c1 <- cima(d$y, d$se, B = 5000, seed = 1, tol = 1e-8, sampling = s)
    p1 <- pima(d$y, d$se, B = 5000, seed = 1, tol = 1e-8, sampling = s)
    expect_identical(c(c1$lci, c1$uci), c(p1$lci, p1$uci))
  }
  c2 <- cima(d$y, d$se, B = 5000, seed = 1, tol = 1e-2, sampling = "exact")
  expect_false(identical(c2$rnd, cima(d$y, d$se, B = 5000, seed = 1,
                                      tol = 1e-8, sampling = "exact")$rnd))
})
