# Prediction intervals with closed-form expressions (methods other than "boot")

pi_formula <- function(muhat, tau2, vmu, nup, nuc, alpha = 0.05) {
  qp <- stats::qt(1 - alpha/2, nup)
  qc <- stats::qt(1 - alpha/2, nuc)
  c(lpi = muhat - qp*sqrt(tau2 + vmu), upi = muhat + qp*sqrt(tau2 + vmu),
    lci = muhat - qc*sqrt(vmu), uci = muhat + qc*sqrt(vmu))
}
lims <- function(r) unlist(r[c("lpi", "upi", "lci", "uci")])

test_that("HTS prediction interval (sbp data; published values)", {
  d <- get_data("sbp")
  res <- pima(d$y, d$se, method = "HTS")
  expect_equal(round(c(res$lpi, res$upi), 4), c(-0.7598, 0.0917))
  expect_equal(round(c(res$lci, res$uci), 4), c(-0.5068, -0.1613))
  expect_equal(round(res$tau2h, 4), 0.0282)
  expect_equal(round(res$i2h, 1), 70.5)
})

test_that("HTS, HK, SJ and APX agree with the formulas", {
  for (d in datasets()) {
    k <- length(d$y)
    v <- d$se^2
    # HTS: DerSimonian-Laird
    t <- ref_tau2(d$y, d$se, "DL")
    w <- 1/(v + t)
    mu <- sum(w*d$y)/sum(w)
    r <- pima(d$y, d$se, method = "HTS")
    expect_equal(r$muhat, mu, tolerance = 1e-10)
    expect_equal(lims(r), pi_formula(mu, t, 1/sum(w), k - 2, k - 1),
                 tolerance = 1e-10, ignore_attr = TRUE)
    # REML based methods (tau^2 from the package to isolate the formulas)
    t <- tau2h(d$y, d$se, method = "REML")$tau2h
    w <- 1/(v + t)
    mu <- sum(w*d$y)/sum(w)
    q <- sum(w*(d$y - mu)^2)/(k - 1)
    h <- w/sum(w)
    vmu <- list(
      HK  = max(1, q)/sum(w),
      SJ  = sum(w^2*(d$y - mu)^2/(1 - h))/sum(w)^2,
      APX = 1/sum(w)
    )
    for (m in names(vmu)) {
      r <- pima(d$y, d$se, method = m)
      expect_equal(r$tau2h, t, tolerance = 1e-12)
      expect_equal(r$vmuhat, vmu[[m]], tolerance = 1e-10)
      expect_equal(lims(r), pi_formula(mu, t, vmu[[m]], k - 2, k - 1),
                   tolerance = 1e-10, ignore_attr = TRUE)
    }
    # "CL" is an alias of "APX"
    expect_equal(lims(pima(d$y, d$se, method = "CL")), lims(pima(d$y, d$se, method = "APX")))
  }
})

test_that("KR: variance inflation, degrees of freedom, and the equal-SE case", {
  for (d in datasets()) {
    r <- pima(d$y, d$se, method = "KR")
    a <- pima(d$y, d$se, method = "APX")
    expect_gte(r$vmuhat, a$vmuhat)
    expect_equal(r$nup, r$nuc - 1)
    expect_gt(r$nuc, 0)
  }
  # with equal standard errors, KR reduces to APX (vmuhat = 1/sum(w), df = K - 1)
  d <- get_data("sbp")
  se <- rep(0.1, length(d$y))
  r <- pima(d$y, se, method = "KR")
  a <- pima(d$y, se, method = "APX")
  expect_equal(r$nuc, length(d$y) - 1, tolerance = 1e-10)
  expect_equal(lims(r), lims(a), tolerance = 1e-10)
})

test_that("intervals: ordering, alpha, and location-scale equivariance", {
  for (m in c("HTS", "HK", "SJ", "KR", "APX")) {
    for (d in datasets()) {
      r <- pima(d$y, d$se, method = m)
      expect_true(r$lpi <= r$lci && r$lci <= r$muhat && r$muhat <= r$uci &&
                    r$uci <= r$upi, label = m)
      r90 <- pima(d$y, d$se, method = m, alpha = 0.1)
      expect_lt(r90$upi - r90$lpi, r$upi - r$lpi)
      s <- pima(3 + 2*d$y, 2*d$se, method = m)
      tol <- if (m == "HTS") 1e-10 else 1e-3  # REML is iterative
      expect_equal((lims(s) - 3)/2, lims(r), tolerance = tol)
    }
  }
})

test_that("HK confidence interval agrees with metafor (test = 'knha')", {
  skip_if_not_installed("metafor")
  for (d in datasets()) {
    fit <- metafor::rma(yi = d$y, sei = d$se, method = "REML", test = "knha")
    r <- pima(d$y, d$se, method = "HK")
    k <- length(d$y)
    q <- sum((d$y - r$muhat)^2/(d$se^2 + r$tau2h))/(k - 1)
    expect_close(r$muhat, as.numeric(fit$beta), tol = 1e-4)
    if (q >= 1.01) {  # pimeta truncates the HK variance at q = 1
      expect_close(r$lci, fit$ci.lb, tol = 1e-3)
      expect_close(r$uci, fit$ci.ub, tol = 1e-3)
    }
  }
})
