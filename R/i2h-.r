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
  
  # initial check
  if (is.null(se)) {
    stop("'se' is a null value.")
  } else if (is.null(tau2h)) {
    stop("'tau2h' is a null value.")
  } else if (any(is.na(se))) {
    stop("'se' has missing value(s).")
  } else if (is.na(tau2h)) {
    stop("'tau2h' is a missing value.")
  } else if (any(is.infinite(se))) {
    stop("'se' has infinite value(s).")
  } else if (is.infinite(tau2h)) {
    stop("'tau2h' is an infinite value.")
  } else if (any(is.nan(se))) {
    stop("'se' has NaN(s).")
  } else if (is.nan(tau2h)) {
    stop("'tau2h' is NaN.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  } else if (tau2h < 0.0) {
    stop("'tau2h' should be positive.")
  }
  
  k <- length(se)
  wi <- se^-2
  s2h <- sum(wi)*(k - 1)/(sum(wi)^2 - sum(wi^2))
  i2h <- 100*tau2h/(s2h + tau2h)
  
  return(i2h)

}