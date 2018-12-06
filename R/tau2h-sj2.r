# Sidik--Jonkman improved estimator for \eqn{\tau^2}
# 
# Returns the Sidik--Jonkman improved estimator
# for \eqn{\tau^2} (Sidik & Jonkman, 2007).
# 
# @name tau2h_sj2
# @rdname tau2h_sj2
# @param y the effect size estimates vector
# @param se the within studies standard errors vector
# @return
# \itemize{
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Sidik, K., and Jonkman, J. N. (2007).
# A comparison of heterogeneity variance estimators
# in combining results of studies.
# \emph{Stat Med.}
# \strong{26}(9): 1964-1981.
# @examples
# data(sbp, package = "pimeta")
# pimeta::tau2h_sj2(sbp$y, sbp$sigmak)
# @export
tau2h_sj2 <- function(y, se) {

  k <- length(y)
  tau2h0 <- max(0.01, tau2h_vc(y, se)$tau2h)
  w <- (1.0 + se^2/tau2h0)^-1
  tau2h <- sum(w*(y - sum(y*w)/sum(w)))/(k - 1.0)
  
  return(list(tau2h = tau2h))
  
}
