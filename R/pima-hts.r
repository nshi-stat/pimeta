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
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' \url{https://doi.org/10.1111/j.1467-985X.2008.00552.x}
#' @seealso
#' \code{\link[=pima]{pima}}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::pima_hts(sbp$y, sbp$sigmak)
#' @export
pima_hts <- function(y, sigma, alpha = 0.05) {

  # initial check
  util_check_num(y)
  util_check_num(sigma)
  util_check_num(alpha)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)

  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  }

  # estimation  
  k <- length(y)
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  vmuhat <- 1/sum(w)
  lpi <- muhat - stats::qt(1 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci,
              tau2h = tau2h, method = "HTS", y = y, se = sigma, alpha = alpha)

  return(res)

}

#' @export
htsdl <- pima_hts

