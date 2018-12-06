# Empirical Bayes estimator for \eqn{\tau^2}
# 
# Returns an Empirical Bayes estimator
# for \eqn{\tau^2} (Morris, 1983).
# 
# @name tau2h_eb
# @rdname tau2h_eb
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @param maxiter the maximum number of iterations
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Morris, C. N. (1983).
# Parametric empirical Bayes inference: theory and applications.
# \emph{J Am Stat Assoc.}
# \strong{78}(381): 47-55.
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_eb(sbp$y, sbp$sigmak)
# @export
tau2h_eb <- function(y, se, maxiter = 100) {

  k <- length(y)
  tau2h <- tau2h_dl(y, se)$tau2h
  r <- 0
  
  # auto ajdustment for step length
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (se^2 + tau2h)^-1
    ti <- sum(wi*(k/(k - 1.0)*(y - sum(wi*y)/sum(wi))^2 - se^2))/sum(wi)
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
  
  return(list(tau2h = tau2h))
  
}
