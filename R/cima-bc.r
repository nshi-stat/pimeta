# Confidence interval based on Bartlett type corrections
#
# Returns a confidence interval for \eqn{\hat{\mu}} based on
# Bartlett type corrections (Noma, 2011).
# 
# @name cima_bc
# @rdname cima_bc
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
# Noma H. (2011)
# Confidence intervals for a random-effects meta-analysis
# based on Bartlett-type corrections.
# \emph{Stat Med.}
# \strong{30}(28): 3304-3312.
# \url{https://doi.org/10.1002/sim.4350}
# @seealso
# \code{\link[=cima]{cima}}.
# @examples
# data(sbp, package = "pimeta")
# pimeta::cima_bc(sbp$y, sbp$sigmak)
# @export
cima_bc <- function(y, se, alpha = 0.05) {

  k <- length(y)
  tau2h <- tau2h_ml(y = y, se = se)$tau2h
  w <- (se^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)

  upper <- muhat + 10*sqrt(sum(w)^-1)
  while(peqnbc(upper, muhat, tau2h, alpha, y, se) > 0) {
    upper <- upper + 10*sqrt(sum(w)^-1)
  }
  res <- try(
    uniroot(peqnbc, interval = c(muhat, upper), muhat = muhat,
            tau2h = tau2h, alpha = alpha, y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  uci <- res$root
  
  lower <- muhat - 10*sqrt(sum(w)^-1)
  while(peqnbc(lower, muhat, tau2h, alpha, y, se) > 0) {
    lower <- lower - 10*sqrt(sum(w)^-1)
  }
  res <- try(
    uniroot(peqnbc, interval = c(lower, muhat), muhat = muhat,
            tau2h = tau2h, alpha = alpha, y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  lci <- res$root
  
  res <- list(muhat = muhat, lci = lci, uci = uci, tau2h = tau2h, vmuhat = NULL,
              method = "BC", y = y, se = se, alpha = alpha)
  
  return(res)
  
}

peqnbc <- function(mu, muhat, tau2h, alpha, y, se) {
  tau2mu <- tau2h_pl(mu, y, se)
  ctau2 <- sum((se^2 + tau2mu)^-3)/(sum((se^2 + tau2mu)^-1)*sum((se^2 + tau2mu)^-2))
  ll(mu, tau2mu, y, se) - ll(muhat, tau2h, y, se) + (0.5 + ctau2)*stats::qchisq(1.0 - alpha, 1)
}

