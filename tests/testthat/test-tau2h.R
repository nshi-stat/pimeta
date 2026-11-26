closed <- c("DL", "VC", "HS", "HM", "SJ", "SJ2")
iterative <- c("PM", "ML", "REML", "AREML", "EB")

test_that("closed-form estimators of tau^2 agree with the formulas", {
  for (d in c(datasets(), list(homog = homog))) {
    for (m in closed) {
      expect_equal(tau2h(d$y, d$se, method = m)$tau2h, ref_tau2(d$y, d$se, m),
                   tolerance = 1e-10, label = m)
    }
  }
})

test_that("iterative estimators of tau^2 solve their estimating equations", {
  for (d in c(datasets(), list(homog = homog))) {
    for (m in iterative) {
      # PM uses uniroot() with the default tolerance (about 1.2e-4)
      tol <- if (m == "PM") 1e-3 else 1e-4
      expect_close(tau2h(d$y, d$se, method = m)$tau2h, ref_tau2(d$y, d$se, m),
                   tol = tol)
    }
  }
})

test_that("estimators of tau^2 are 0 for homogeneous data", {
  for (m in c(closed, iterative)) {
    if (m %in% c("SJ", "SJ2", "HM")) next  # these are positive when Q > 0
    expect_equal(tau2h(homog$y, homog$se, method = m)$tau2h, 0, label = m)
  }
  expect_gt(tau2h(homog$y, homog$se, method = "SJ")$tau2h, 0)
  expect_gt(tau2h(homog$y, homog$se, method = "BM")$tau2h, 0)
  expect_gt(tau2h(homog$y, homog$se, method = "HM")$tau2h, 0)
})

test_that("estimators of tau^2 are non-negative and the BM estimate exceeds ML", {
  for (d in datasets()) {
    for (m in c(closed, iterative, "BM")) {
      expect_gte(tau2h(d$y, d$se, method = m)$tau2h, 0)
    }
    ml <- tau2h(d$y, d$se, method = "ML")$tau2h
    expect_gt(tau2h(d$y, d$se, method = "BM")$tau2h, ml)
  }
})

test_that("estimators of tau^2 are scale equivariant", {
  d <- get_data("sbp")
  for (m in c(closed, iterative, "BM")) {
    if (m == "SJ2" || m == "BM") next  # these use absolute constants
    a <- tau2h(d$y, d$se, method = m)$tau2h
    b <- tau2h(10*d$y, 10*d$se, method = m)$tau2h
    expect_close(b/100, a, tol = if (m %in% iterative) 1e-3 else 1e-10)
  }
})

test_that("Wald confidence intervals for tau^2", {
  z <- stats::qnorm(0.975)
  for (d in datasets()) {
    v <- d$se^2
    t <- tau2h(d$y, d$se, method = "ML")$tau2h
    w <- 1/(v + t)
    se <- sqrt(2/sum(w^2))
    r <- tau2h(d$y, d$se, method = "ML", methodci = "ML")
    expect_equal(c(r$lci, r$uci), c(max(0, t - z*se), t + z*se), tolerance = 1e-10)
    t <- tau2h(d$y, d$se, method = "REML")$tau2h
    w <- 1/(v + t)
    # REML information: tr(P^2)/2 with P = W - W 1 1' W / sum(w)
    P <- diag(w) - outer(w, w)/sum(w)
    se <- sqrt(2/sum(P*P))
    r <- tau2h(d$y, d$se, method = "DL", methodci = "REML")
    expect_equal(c(r$lci, r$uci), c(max(0, t - z*se), t + z*se), tolerance = 1e-10)
    # the point estimate is the one of 'method'
    expect_equal(r$tau2h, ref_tau2(d$y, d$se, "DL"), tolerance = 1e-10)
  }
})

test_that("tau^2 estimators agree with metafor", {
  skip_if_not_installed("metafor")
  mf <- c(DL = "DL", VC = "HE", HS = "HS", SJ = "SJ", PM = "PM",
          ML = "ML", REML = "REML", EB = "EB")
  ds <- datasets()
  for (nm in names(ds)) {
    d <- ds[[nm]]
    for (m in names(mf)) {
      est <- tau2h(d$y, d$se, method = m)$tau2h
      fit <- metafor::rma(yi = d$y, sei = d$se, method = mf[[m]])
      if (nm == "hyp" && m == "ML") {
        # The ML log-likelihood of 'hyp' has a local maximum at tau^2 = 0
        # (the score is negative at 0) besides the global maximum at about
        # 1.076; metafor's Fisher scoring stops at the boundary. Check that
        # the pimeta estimate is the global maximiser instead.
        expect_close(est, ref_tau2(d$y, d$se, "ML"), tol = 1e-3)
        lml <- function(t) {
          w <- 1/(d$se^2 + t)
          ref_ll(sum(w*d$y)/sum(w), t, d$y, d$se^2)
        }
        expect_gt(lml(est), lml(fit$tau2))
        expect_close(est, metafor::rma(yi = d$y, sei = d$se, method = "ML",
                                       control = list(tau2.init = 1))$tau2, tol = 1e-3)
        next
      }
      tol <- if (m %in% iterative) 1e-3 else 1e-8
      expect_close(est, fit$tau2, tol = tol)
    }
    fit <- metafor::rma(yi = d$y, sei = d$se, method = "REML")
    r <- tau2h(d$y, d$se, method = "REML", methodci = "REML")
    z <- stats::qnorm(0.975)
    expect_close(r$uci, fit$tau2 + z*fit$se.tau2, tol = 1e-3)
  }
})

test_that("print.pima_tau2h", {
  d <- get_data("sbp")
  for (m in c(closed, iterative, "BM")) {
    expect_output(print(tau2h(d$y, d$se, method = m)), "Heterogeneity Variance")
  }
  expect_output(print(tau2h(d$y, d$se, methodci = "REML")), "confidence interval")
})
