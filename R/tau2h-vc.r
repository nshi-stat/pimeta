#' Variance component type estimator for \eqn{\tau^2}
#' 
#' Returns the variance component type estimator
#' for \eqn{\tau^2} (Hedges, 1983).
#' 
#' @name tau2h_vc
#' @rdname tau2h_vc
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Hedges, L. V. (1983).
#' A random effects model for effect sizes.
#' \emph{Psychol Bull.}
#' \strong{93}(2): 388-395. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h_vc(sbp$y, sbp$sigmak)
#' @export
tau2h_vc <- function(y, se) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  k <- length(y)
  tau2h <- max(0.0, sum(y - sum(y)/k)^2/(k - 1.0) - sum(se)/k)
  
  return(tau2h)
  
}
