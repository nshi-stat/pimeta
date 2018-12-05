#' Paule--Mandel estimator for \eqn{\tau^2}
#' 
#' A subroutine of the Paule--Mandel estimator
#' for \eqn{\tau^2} (Paule & Mandel, 1982).
#' 
#' @name tau2h_pm
#' @rdname tau2h_pm
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Mandel, R. C., and Mandel, K. H. (1982).
#' Consensus values and weighting factors. 
#' \emph{J Res Natl Inst Stand Techno.}
#' \strong{87}(5): 377-385. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h_pm(sbp$y, sbp$sigmak)
#' @export
tau2h_pm <- function(y, se) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  qgen <- function(x, y, se) {
    k <- length(y)
    wi <- (se^2 + x)^-1
    muhat <- sum(wi*y)/sum(wi)
    qgen <- sum(wi * (y - muhat)^2) - (k - 1)
  }
  
  if (qgen(0, yi, se) < 0) {
    tau2h <- 0  
  } else {
    upper <- 20
    while(qgen(upper, yi, se) > 0) {
      upper <- upper*2
    }
    res <- try(
      uniroot(qgen, interval = c(0, upper), y = yi, se = sqrt(vi)),
      silent = FALSE
    )
    if (class(res) == "try-error") {
      stop("Error in tau2h_pm() :\nTry another estimator.")
      res <- NA
    }
    tau2h <- res$root
  }

  return(tau2h)
  
}