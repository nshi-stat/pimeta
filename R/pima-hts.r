#' @export
pima_hts <- function(y, sigma, alpha = 0.05) {

  # initial check
  util_check_num(y)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)

  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  }

  # estimation
  k <- length(y)
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  vmuhat <- 1/sum(w)
  lpi <- muhat - stats::qt(1 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci,
              tau2h = tau2h, vmuhat = vmuhat, nup = k - 2, nuc = k - 1,
              method = "HTS", y = y, se = sigma, alpha = alpha)

  return(res)

}

#' @export
htsdl <- pima_hts

