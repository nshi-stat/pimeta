#' @export
pima_boot <- function(y, sigma, alpha = 0.05, B = 25000, maxit1 = 100000,
                      eps = 10^(-10), lower = 0, upper = 1000, maxit2 = 1000,
                      tol = .Machine$double.eps^0.25, rnd = NULL, parallel = FALSE,
                      seed = NULL) {
  
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
    set.seed(seed)
    rndtau2 <- rtau2CppWrap(
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
      nthread = as.integer(parallel)
    )
  } else {
    rndtau2 <- rnd
  }
  
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- list(muhat = sum(y*w) / sum(w))
  res <- bootPICppWrap(
    rnd   = as.vector(rndtau2),
    y     = as.vector(y),
    sigma = as.vector(sigma),
    alpha = as.double(alpha)
  )
  
  res <- append(append(muhat, res),
                list(tau2h = tau2h, vmuhat = NULL, nup = k - 1, nuc = k - 1,
                     method = "boot", y = y, se = sigma,
                     alpha = alpha, rnd = rndtau2))
  
  return(res)
  
}


#' @export
bootPI <- pima_boot
