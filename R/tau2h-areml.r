# Approximate restricted maximum likelihood estimator for \eqn{\tau^2}
# 
# Returns an approximate restricted maximum likelihood estimator
# for \eqn{\tau^2} (e.g., Thompson & Sharp, 1999).
# 
# @name tau2h_areml
# @rdname tau2h_areml
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @param maxiter the maximum number of iterations
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Thompson, S. G., and Sharp, S. J. (1999).
# Explaining heterogeneity in meta-analysis: a comparison of methods.
# \emph{Stat Med.}
# \strong{18}(20): 2693-2708.
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_areml(sbp$y, sbp$sigmak)
# @export
tau2h_areml <- function(y, se, maxiter = 100) {

  tau2h <- tau2h_dl(y, se)
  r <- 0
  k <- length(y)
  
  # auto ajdustment for step length
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (se^2 + tau2h)^-1
    ti <- sum(wi^2*((k/(k - 1))*(y - sum(wi*y)/sum(wi))^2 - se^2))/sum(wi^2)
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
