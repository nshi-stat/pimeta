#' @export
pima_htsreml <- function(y, sigma, alpha = 0.05,
                         vartype = c("HK", "SJBC", "KR", "CL", "APX"), maxiter = 100) {

  # initial check
  lstm <- c("HK", "SJBC", "KR", "CL", "APX")
  vartype <- match.arg(vartype)
  
  util_check_num(y)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (!is.element(vartype, lstm)) {
    stop("Unknown 'vartype' specified.")
  }
  
  # estimation
  k <- length(y)
  tau2h <- tau2h_reml(y = y, se = sigma, maxiter = maxiter)$tau2h
  
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  w1p <- sum(w)
  w2p <- sum(w^2)
  w3p <- sum(w^3)
  nup <- k - 2
  nuc <- k - 1
  
  if (vartype == "HK") {
    q <- sum(w*(y - muhat)^2)/(k - 1.0)
    vmuhat <- max(1.0, q)/w1p
    method <- "HK"
  } else if (vartype == "SJBC") {
    lambda <- sigma^2 + tau2h
    h <- 2.0*w/w1p - sum(w^2*lambda)/(lambda*w1p^2)
    vmuhat <- sum(w^2*(y - muhat)^2/(1.0 - h))/w1p^2
    method <- "SJ"
  } else if (vartype == "KR") {
    IE <- w2p*0.5 - w3p/w1p + 0.5*(w2p/w1p)^2
    vmuhat <- (1.0 + 2.0*(w3p/w1p - (w2p/w1p)^2)/IE)/w1p
    nup <- 2*IE/(vmuhat*w2p)^2 - 1
    nuc <- 2*IE/(vmuhat*w2p)^2
    method <- "KR"
  } else if (vartype == "CL" | vartype == "APX") {
    vmuhat <- 1.0/w1p
    method <- "APX"
  } else {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/w1p
    method <- "HK"
  }
  
  if (nup <= 0) {
    lpi <- upi <- lci <- uci <- NaN
  } else {
    lpi <- muhat - stats::qt(1.0 - alpha*0.5, nup)*sqrt(tau2h + vmuhat)
    upi <- muhat + stats::qt(1.0 - alpha*0.5, nup)*sqrt(tau2h + vmuhat)
    lci <- muhat - stats::qt(1.0 - alpha*0.5, nuc)*sqrt(vmuhat)
    uci <- muhat + stats::qt(1.0 - alpha*0.5, nuc)*sqrt(vmuhat)
  }
  
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci, 
              tau2h = tau2h, vmuhat = vmuhat, nup = nup, nuc = nuc,
              method = method, y = y, se = sigma, alpha = alpha)

  return(res)

}


#' @export
htsreml <- pima_htsreml
