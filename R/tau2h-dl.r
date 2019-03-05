# DerSimonian--Laird estimator for \eqn{\tau^2}
# 
# Returns the DerSimonian--Laird estimator
# for \eqn{\tau^2} (DerSimonian & Laird, 1986).
# 
# @name tau2h_dl
# @rdname tau2h_dl
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
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
# pimeta::tau2h_dl(sbp$y, sbp$sigmak)
# @export
tau2h_dl <- function(y, se) {

  k <- length(y)
  vi <- se^-2
  Q <- sum(vi*(y - sum(vi*y)/sum(vi))^2)
  tau2h <- max(0.0, (Q - (k - 1.0)) /
                 (sum(vi) - sum(vi^2)/sum(vi)))
  
  return(list(tau2h = tau2h, Q = Q))
  
}
