#' Higgins-Thompson-Spiegelhalter prediction interval with the Dersimonian-Laird estimator for \eqn{\hat{\tau}}
#'
#' A subroutine for the original HTS PI
#' 
#' @name pima_hts
#' @rdname pima_hts
#' @aliases htsdl
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @return The average treatment effect estimate \eqn{\hat{\mu}} (\code{muhat}),
#' the prediction limits \eqn{\hat{c}_l} (\code{lpi}) and \eqn{\hat{c}_u} (\code{upi}),
#' the confidence limits \eqn{\hat{\mu}_l} (\code{lci}) and \eqn{\hat{\mu}_u} (\code{uci}),
#' the DL estimator for \eqn{\hat{\tau}}
#' (\code{tau2h}).
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::pima_hts(sbp$y, sbp$sigmak)
#' # $muhat
#' # [1] -0.3340597
#' # $lbpi
#' # [1] -0.7597777
#' # $ubpi
#' # [1] 0.09165839
#' # $tau2h
#' # [1] 0.02824971
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
  res <- list(muhat = muhat, lpi = lpi, upi = upi,
              lci = lci, uci = uci, tau2h = tau2h)

  return(res)

}

#' @export
htsdl <- pima_hts

