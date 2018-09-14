#' Higgins--Thompson--Spiegelhalter prediction interval
#'
#' A subroutine for the Higgins--Thompson--Spiegelhalter PI
#' based on the DerSimonian-Laird estimator (Higgins et al., 2009)
#' 
#' @name pima_hts
#' @rdname pima_hts
#' @aliases htsdl
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' @seealso
#' \code{\link[=pima]{pima()}}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::pima_hts(sbp$y, sbp$sigmak)
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # Higgins-Thompson-Spiegelhalter prediction interval
#' #  Heterogeneity variance: DerSimonian-Laird
#' #  SE for average treatment effect: standard
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3341 [-0.7598, 0.0917]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3341 [-0.5068, -0.1613]
#' # 
#' # Heterogeneity variance (tau^2):
#' #  0.0282
#' # 
#' @export
pima_hts <- function(y, sigma, alpha = 0.05) {

  ## .. need more more strictry check.
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (min(sigma) < 0.0) {
    stop("'sigma' should be positive.")
  }

  k <- length(y)
  tau2h <- max(0, (sum(sigma^-2 * (y - sum(sigma^-2*y) / sum(sigma^-2))^2) - (k - 1)) /
                 (sum(sigma^-2) - sum(sigma^-4)/sum(sigma^-2)))
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  vmuhat <- 1/sum(w)
  lpi <- muhat - stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci,
              tau2 = tau2h, method = "HTS", y = y, se = sigma, alpha = alpha)
  class(res) <- "pima" 
  
  return(res)

}

#' @export
htsdl <- pima_hts

