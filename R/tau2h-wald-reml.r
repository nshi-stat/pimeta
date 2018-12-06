# Wald REML confidence interval for \eqn{\tau^2}
# 
# Returns a Wald confidence interval for \eqn{\tau^2}
# with REML estimator (Veroniki et al., 2016).
# 
# @name tau2h_wald_reml
# @rdname tau2h_wald_reml
# @param se the within studies standard errors vector
# @param tau2h the estimate of \eqn{\tau^2}
# @param alpha the alpha level of the confidence interval
# @return
# \itemize{
# \item \code{lci}, \code{uci}: the lower and upper confidence limits
#       \eqn{\hat{\tau}^2_l} and \eqn{\hat{\tau}^2_u}.
# }
# @references
# Veroniki, A. A., Jackson, D., Viechtbauer, W.,
# Bender, R., Bowden, J., Knapp, G., Kuss, O.,
# Higgins, J. P. T., Langan, D., and Salanti, J. (2016).
# Methods to estimate the between-study variance
# and its uncertainty in meta-analysis.
# \emph{Res Syn Meth.}
# \strong{7}(1): 55-79.
# @examples
# data(sbp, package = "pimeta")
# tau2h <- pimeta::tau2h_reml(sbp$y, sbp$sigmak)
# pimeta::tau2h_wald_reml(sbp$sigmak, tau2h)
# @export
tau2h_wald_reml <- function(se, tau2h, alpha = 0.05) {
  
  # initial check
  if (is.null(tau2h)) {
    stop("'tau2h' is a null value.")
  } else if (is.na(tau2h)) {
    stop("'tau2h' is a missing value.")
  } else if (is.infinite(tau2h)) {
    stop("'tau2h' is an infinite value.")
  } else if (is.nan(tau2h)) {
    stop("'tau2h' is NaN.")
  } else if (tau2h < 0.0) {
    stop("'tau2h' should be positive.")
  }
  
  wi <- (se^2 + tau2h)^-1
  vtau2h <- 2.0*(sum(wi^2) - 2.0*sum(wi^3)/sum(wi) + sum(wi^2)^2/sum(wi)^2)^-1
  lci <- muhat - stats::qnorm(1.0 - alpha*0.5)*sqrt(vtau2h)
  uci <- muhat + stats::qnorm(1.0 - alpha*0.5)*sqrt(vtau2h)
  res <- list(lci = lci, uci = uci, alpha = alpha)
  
  return(res)

}
