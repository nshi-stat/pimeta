#' @export
pima_boot <- function(y, sigma, alpha = 0.05, B = 25000, maxit1 = 100000,
                      eps = 10^(-10), lower = 0, upper = 1000, maxit2 = 1000,
                      tol = .Machine$double.eps^0.25, rnd = NULL, parallel = FALSE,
                      seed = NULL, sampling = c("table", "exact"), nodes = 400) {
  
  sampling <- match.arg(sampling)
  
  # initial check
  util_check_num(y)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  util_check_gt(B, 1)
  util_check_nonneg(parallel)
  util_check_gt(maxit1, 1)
  util_check_gt(eps, 0)
  util_check_ge(lower, 0)
  util_check_gt(upper, 0)
  util_check_gt(maxit2, 1)
  util_check_gt(tol, 0)
  util_check_gt(nodes, 3)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (lower >= upper) {
    stop("'upper' should be greater than 'lower'.")
  }
  
  if (B < 5000) {
    warning("'B' > 5000 is recommended.")
  }
  
  k <- length(y)
  
  # random numbers generation
  if (is.null(rnd)) {
    if (parallel == FALSE) {
      parallel <- 1
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }
    rndtau2 <- rtau2CDCppWrap(
      n       = as.integer(B),
      y       = as.vector(y),
      sigma   = as.vector(sigma),
      mode    = as.double(1),
      maxit1  = as.integer(maxit1),
      eps     = as.double(eps),
      lower   = as.double(lower),
      upper   = as.double(upper),
      maxit2  = as.integer(maxit2),
      tol     = as.double(tol),
      nthread = as.integer(parallel),
      table   = as.integer(sampling == "table"),
      nodes   = as.integer(nodes)
    )
  } else {
    rndtau2 <- rnd
  }
  
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- list(muhat = sum(y*w) / sum(w))
  
  # tau^2 = Inf (beyond the upper end of the search) is handled as the limit
  # tau^2 -> Inf: equal weights and a predictive draw of +-Inf. It is passed to
  # bootPICppWrap() as a huge finite value (with Inf, the weights are 0 and the
  # draws become NaN, which made the quantiles unreliable). Draws of tau^2 that
  # could not be generated (NA; see the warning of rtau2CDCppWrap) are dropped.
  ninf <- sum(is.infinite(rndtau2))
  rnduse <- as.vector(rndtau2)[!is.na(rndtau2)]
  rnduse[is.infinite(rnduse)] <- 1e200
  res <- bootPICppWrap(
    rnd   = rnduse,
    y     = as.vector(y),
    sigma = as.vector(sigma),
    alpha = as.double(alpha)
  )
  if (ninf > 0) {
    toinf <- function(x) {
      big <- is.finite(x) & abs(x) >= 1e90
      x[big] <- sign(x[big])*Inf
      x
    }
    res$lpi <- toinf(res$lpi)
    res$upi <- toinf(res$upi)
    warning(sprintf("%d random numbers of tau^2 exceeded the upper limit of the search and were treated as tau^2 = Inf.", ninf))
  }
  if (is.null(rnd) && length(unique(y)) == 1) {
    warning("All the estimates 'y' are identical: all random numbers of tau^2 are 0 and the intervals have zero width.")
  }
  
  res <- append(append(muhat, res),
                list(tau2h = tau2h, vmuhat = NULL, nup = k - 1, nuc = k - 1,
                     method = "boot", y = y, se = sigma,
                     alpha = alpha, rnd = rndtau2,
                     sampling = if (is.null(rnd)) sampling else "user"))
  
  return(res)
  
}


#' @export
bootPI <- pima_boot
