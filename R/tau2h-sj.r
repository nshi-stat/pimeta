#' Sidik--Jonkman estimator for \eqn{\tau^2}
#' 
#' A subroutine of the Sidik--Jonkman estimator
#' for \eqn{\tau^2} (Sidik & Jonkman, 2005).
#' 
#' @name tau2h_sj
#' @rdname tau2h_sj
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Sidik, K., and Jonkman, J. N. (2005).
#' Simple heterogeneity variance estimation for meta-analysis.
#' \emph{J R Stat Soc Ser C Appl Stat.}
#' \strong{54}(2): 367-384. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h_sj(sbp$y, sbp$sigmak)
#' @export
tau2h_sj <- function(y, se) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  k <- length(y)
  tau2h0 <- sum((y - sum(y)/k)^2)/k
  w <- (1.0 + se^2/tau2h0)^-1
  tau2h <- max(0.0, sum(w*(y - sum(y*w)/sum(w)))/(k - 1.0))
  
  return(tau2h)
  
}
