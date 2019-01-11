# Profile likelihood confidence interval
#
# Returns a profile likelihood confidence interval
# for \eqn{\hat{\mu}} (Hardy & Thompson, 1996).
# 
# @name cima_pl
# @rdname cima_pl
# @param y the effect size estimates vector
# @param sigma the within studies standard errors vector
# @param alpha the alpha level of the prediction interval
# @return
# \itemize{
# \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
# \item \code{lci}, \code{lci}: the lower and upper confidence limits
#                               \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Hardy, R. J., and Thompson, S. G. (1996).
# A likelihood approach to meta-analysis with random effects.
# \emph{Stat Med.}
# \strong{15}(6): 619-629.
# \url{https://doi.org/10.1002/(SICI)1097-0258(19960330)15:6<619::AID-SIM188>3.0.CO;2-A}.
# @seealso
# \code{\link[=cima]{cima}}.
# @examples
# data(sbp, package = "pimeta")
# pimeta::cima_pl(sbp$y, sbp$sigmak)
# @export
cima_pl <- function(y, se, alpha = 0.05) {

  k <- length(y)
  tau2h <- tau2h_ml(y = y, se = se)$tau2h
  w <- (se^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  
  upper <- muhat + 10*sqrt(sum(w)^-1)
  while(peqn(upper, muhat, tau2h, alpha, y, se) > 0) {
    upper <- upper + 10*sqrt(sum(w)^-1)
  }
  res <- try(
    uniroot(peqn, interval = c(muhat, upper), muhat = muhat,
            tau2h = tau2h, alpha = alpha, y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  uci <- res$root
  
  lower <- muhat - 10*sqrt(sum(w)^-1)
  while(peqn(lower, muhat, tau2h, alpha, y, se) > 0) {
    lower <- lower - 10*sqrt(sum(w)^-1)
  }
  res <- try(
    uniroot(peqn, interval = c(lower, muhat), muhat = muhat,
            tau2h = tau2h, alpha = alpha, y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  lci <- res$root
  
  res <- list(muhat = muhat, lci = lci, uci = uci, tau2h = tau2h,
              method = "PL", y = y, se = se, alpha = alpha)
  
  return(res)
  
}

peqn <- function(mu, muhat, tau2h, alpha, y, se) {
  tau2mu <- tau2h_pl(mu, y, se)
  ll(mu, tau2mu, y, se) - ll(muhat, tau2h, y, se) + 0.5*stats::qchisq(1.0 - alpha, 1)
}

ll <- function(mu, tau2, y, se) {
  -0.5*sum(log(se^2 + tau2)) - 0.5*sum((y - mu)^2/(se^2 + tau2))
}

tau2h_pl <- function(mu, y, se, maxiter = 100) {
  
  tau2h <- tau2h_dl(y, se)$tau2h
  r <- 0
  
  # auto ajdustment for step length
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (se^2 + tau2h)^-1
    ti <- sum(wi^2 * ((y - mu)^2 - se^2)) / sum(wi^2)
    if (ti <= 0) {
      tau2h <- 0.0
      break
    } else {
      if (abs(ti - tau2h)/(1.0 + tau2h) < 1e-5) {
        # converged
        break
      } else if (r == maxiter && autoadj == 1) {
        # not converged
        stop("The heterogeneity variance (tau^2) could not be calculated.")
        break
      } else if (r == maxiter && autoadj == 0) {
        # not converged but retry with an adjusted step length
        r <- 0
        autoadj <- 1
      } else if (autoadj == 1) {
        # with the adjustment
        r <- r + 1
        tau2h <- tau2h + (ti - tau2h)*stepadj
      } else {
        # without the adjustment
        r <- r + 1
        tau2h <- ti
      }
    }
  }
  
  return(tau2h)

}

