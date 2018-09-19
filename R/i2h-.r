#' \eqn{I^2} heterogeneity measure
#' 
#' A subroutine of a estimator for
#' \eqn{I^2} (Higgins & Thompson, 2002).
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
#' @examples
#' data(sbp, package = "pimeta")
#' tau2h <- pimeta::tau2h_dl(sbp$y, sbp$sigmak)
#' pimeta::i2h(sbp$sigmak, tau2h)
#' @export
i2h <- function(se, tau2h) {
  
  k <- length(se)
  wi <- se^-2
  s2h <- sum(wi)*(k - 1)/(sum(wi)^2 - sum(wi^2))
  i2h <- 100*tau2h/(s2h + tau2h)
  
  return(i2h)

}