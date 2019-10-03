pima_cprob <- function(x, theta0, side) {
  
  if (x$method == "boot") {
    if (side == "lt") {
      cprob <- mean(x$bspi < theta0)
    } else {
      cprob <- mean(x$bspi > theta0)
    }
  } else {
    f <- function(x, nup, muhat, vmuhat, tau2h) {
      theta0 + muhat + qt(x, nup)*sqrt(vmuhat + tau2h)
    }
    cprob <- NULL
    if (x$muhat >= theta0) {
      cprob <- uniroot(f, c(0, 0.5), nup = x$nup, muhat = x$muhat,
                       vmuhat = x$vmuhat, tau2h = x$tau2h)
    } else if (x$muhat < 0) {
      cprob <- uniroot(f, c(0.5, 1), nup = x$nup, muhat = x$muhat,
                       vmuhat = x$vmuhat, tau2h = x$tau2h)
    }
    if (side == "lt") {
      cprob <- ifelse(is.null(cprob), 0, cprob$root)
    } else {
      cprob <- ifelse(is.null(cprob), 1, 1 - cprob$root)
    }
  }
  return(cprob)
  
}
