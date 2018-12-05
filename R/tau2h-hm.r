#' Hartung--Makambi estimator for \eqn{\tau^2}
#' 
#' Returns the Hartung--Makambi estimator
#' for \eqn{\tau^2} (Hartung & Makambi, 1983).
#' 
#' @name tau2h_hm
#' @rdname tau2h_hm
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Hartung, J., and Makambi, K. H. (2003).
#' Reducing the number of unjustified significant results
#' in meta-analysis. 
#' \emph{Commun Stat Simul Comput.}
#' \strong{32}(4): 1179-1190. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h_hm(sbp$y, sbp$sigmak)
#' @export
tau2h_hm <- function(y, se) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  k <- length(y)
  qhm <- sum(se^-2*(y - sum(y*se^-2)/sum(se^-2))^2)
  tau2h <- max(0.0, qhm^2/((sum(se^-2) - sum(se^-4)/sum(se^-2))*(2.0*(k - 1.0) + qhm)))
  
  return(tau2h)
  
}
