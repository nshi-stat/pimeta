# Maximum likelihood estimator for \eqn{\tau^2}
# 
# Returns a maximum likelihood estimator
# for \eqn{\tau^2} (e.g., DerSimonian & Laird, 1986).
# 
# @name tau2h_ml
# @rdname tau2h_ml
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @param maxiter the maximum number of iterations
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# DerSimonian, R., and Laird, N. (1986).
# Meta-analysis in clinical trials.
# \emph{Control Clin Trials.}
# \strong{7}(3): 177-188. 
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_ml(sbp$y, sbp$sigmak)
# @export
tau2h_ml <- function(y, se, maxiter = 100) {

  tau2h <- tau2h_dl(y, se)$tau2h
  r <- 0
  
  # auto ajdustment for step length
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (se^2 + tau2h)^-1
    ti <- sum(wi^2 * ((y - sum(wi*y)/sum(wi))^2 - se^2)) / sum(wi^2)
    if (ti <= 0) {
      tau2h <- 0.0
      break
    } else {
      if (abs(ti - tau2h)/(1.0 + tau2h) < 1e-5) {
        # converged
        break
      } else if (r == maxiter && autoadj == 1) {
        # not converged
        warning("The heterogeneity variance (tau^2) could not be calculated.")
        tau2h <- NaN
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
