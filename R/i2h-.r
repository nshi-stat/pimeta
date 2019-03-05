#' \eqn{I^2} heterogeneity measure
#' 
#' Returns the estimator for
#' (Higgins & Thompson, 2002).
#' 
#' @name i2h
#' @rdname i2h
#' @param se the within studies standard errors vector
#' @param tau2h the estimate of \eqn{\tau^2}
#' @return
#' \itemize{
#' \item \code{i2h}: the estimate for \eqn{I^2}.
#' }
#' @references
#' Higgins, J. P. T., and Thompson, S. G. (2002).
#' Quantifying heterogeneity in a meta-analysis.
#' \emph{Stat Med.}
#' \strong{21}(11): 1539-1558.
#' \url{https://doi.org/10.1002/sim.1186}
#' @examples
#' data(sbp, package = "pimeta")
#' tau2h <- pimeta::tau2h(sbp$y, sbp$sigmak)
#' pimeta::i2h(sbp$sigmak, tau2h$tau2h)
#' @export
i2h <- function(se, tau2h) {
  
  # initial check
  util_check_num(se)
  util_check_num(tau2h)
  util_check_nonneg(se)
  
  # estimation
  k <- length(se)
  wi <- se^-2
  s2h <- sum(wi)*(k - 1)/(sum(wi)^2 - sum(wi^2))
  i2h <- 100*tau2h/(s2h + tau2h)
  
  return(i2h)

}