# Bayes modal estimator for \eqn{\tau^2}
# 
# Returns the Chung--Rabe-Hesketh--Choi Bayes modal
# estimator with a gamma prior \eqn{G(2, 10^{-4})}
# for \eqn{\tau^2} (Chung et al., 2013).
# 
# @name tau2h_bm
# @rdname tau2h_bm
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @param maxiter the maximum number of iterations
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Chung, Y. L., Rabe-Hesketh, S., and Choi, I-H. (2013).
# Avoiding zero between-study variance estimates
# in random-effects meta-analysis.
# \emph{Stat Med.}
# \strong{32}(23): 4071-4089.
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_bm(sbp$y, sbp$sigmak)
# @export
tau2h_bm <- function(y, se, maxiter = 100) {

  tau2h <- tau2h_ml(y, se, maxiter = maxiter)
  if (tau2h > 0) {
    wi <- (se^2 + tau2h)^-1
    vml <- 2.0*sum(wi^2)^-1
    tau2h <- (sqrt(tau2h)*0.5 + sqrt(tau2h)*0.5*sqrt(1 + 4*vml/tau2h))^2
  } else {
    vi <- se^-2
    vml <- 2.0*sum(vi^2)^-1
    tau2h <- vml
  }

  return(tau2h)
  
}
