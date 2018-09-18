#' DerSimonian--Laird estimator for \eqn{\tau^2}
#' 
#' A subroutine of the DerSimonian--Laird estimator
#' for \eqn{\tau^2} (DerSimonian & Laird, 1986).
#' 
#' @name tau2h_dl
#' @rdname tau2h_dl
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' DerSimonian, R., and Laird, N. (1986).
#' Meta-analysis in clinical trials.
#' \emph{Control Clin Trials.}
#' \strong{7}(3): 177-188. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h_dl(sbp$y, sbp$sigmak)
#' @export
tau2h_dl <- function(y, se) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  k <- length(y)
  tau2h <- max(0.0, (sum(se^-2*(y - sum(se^-2*y)/sum(se^-2))^2) - (k - 1.0)) /
                 (sum(se^-2) - sum(se^-4)/sum(se^-2)))
  
  return(tau2h)
  
}
