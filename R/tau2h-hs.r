# Hunter--Schmidt estimator for \eqn{\tau^2}
# 
# Returns the Hunter--Schmidt estimator
# for \eqn{\tau^2} (Hunter & Schmidt, 2004).
# This estimator has negative bias (Viechtbauer, 2005).
# 
# @name tau2h_hs
# @rdname tau2h_hs
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Hunter, J. E., and Schmidt, F. L. (2004).
# \emph{Methods of Meta-Analysis: Correcting Error and Bias in Research Findings. 2nd edition.}
# Sage Publications, Inc.
# 
# Viechtbauer, W. (2005).
# Bias and efficiency of meta-analytic variance
# estimators in the random-effects Model.
# \emph{J Educ Behav Stat.}
# \strong{30}(3): 261-293.
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_hs(sbp$y, sbp$sigmak)
# @export
tau2h_hs <- function(y, se) {

  k <- length(y)
  vi <- se^-2
  tau2h <- max(0.0, (sum(vi*(y - sum(vi*y)/sum(vi))^2) - (k - 1.0)) / (sum(vi)))
  
  return(tau2h)
  
}
