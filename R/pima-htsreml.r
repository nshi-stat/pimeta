#' Higgins-Thompson-Spiegelhalter prediction interval with a REML estimator for \eqn{\hat{\tau}}
#'
#' A subroutine for the extended HTS PI (Partlett & Riley, 2017)
#' 
#' @name pima_htsreml
#' @rdname pima_htsreml
#' @aliases htsreml
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @param vartype the type of the variance estimator for \eqn{\hat{\mu}} (default = "HK"):
#' \code{HK}, the Hartung and Knapp (2001)'s estimator;
#' \code{SJBC}, the Sidik and Jonkman (2006)'s bias coreccted estimator;
#' \code{CL}, a classical estimator, \eqn{(1/\sum{\hat{w}_i})^{-1}};
#' @param maxiter the maximum number of iterations
#' @return The average treatment effect estimate \eqn{\hat{\mu}} (\code{muhat}),
#' the lower and upper prediction limits \eqn{\hat{c}_l} (\code{lpi}) and \eqn{\hat{c}_u} (\code{upi}),
#' the REML estimator for \eqn{\hat{\tau}} (\code{tau2h}), and the type of the variance estimator (\code{vartype}).
#' @references
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::htsreml(sbp$y, sbp$sigmak)
#' # $muhat
#' # [1] -0.3287403
#' # $lbpi
#' # [1] -0.9886995
#' # $ubpi
#' # [1] 0.3312188
#' # $tau2h
#' # [1] 0.06995102
#' @export
pima_htsreml <- function(y, sigma, alpha = 0.05,
                         vartype = c("HK", "SJBC", "CL"), maxiter = 100) {

  vartype <- match.arg(vartype)

  ## .. need more more strictry check.
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (min(sigma) < 0.0) {
    stop("'sigma' should be positive.")
  }

  k <- length(y)
  taudl <- max(0, (sum(sigma^-2 * (y - sum(sigma^-2*y) / sum(sigma^-2))^2) - (k - 1)) /
                 (sum(sigma^-2) - sum(sigma^-4)/sum(sigma^-2)))
  tau2h <- taudl
  r <- 0
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (sigma^2 + tau2h)^-1
    ti <- sum(wi^2 * ((y - sum(wi*y)/sum(wi))^2 + 1/sum(wi) - sigma^2)) / sum(wi^2)
    if (ti <= 0) {
      tau2h <- 0
      break
    } else {
      if (abs(ti - tau2h)/(1 + tau2h) < 1e-5) {
        break
      } else if (r == maxiter && autoadj == 1) {
        break
      } else if (r == maxiter && autoadj == 0) {
        r <- 0
        autoadj <- 1
      } else if (autoadj == 1) {
        r <- r + 1
        tau2h <- tau2h + (ti - tau2h)*stepadj
      } else {
        r <- r + 1
        tau2h <- ti
      }
    }
  }
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)

  if (vartype == "HK") {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/sum(w)
  } else if (vartype == "SJBC") {
    lambda <- sigma^2 + tau2h
    h <- 2*w/sum(w) - sum(w^2*lambda)/lambda/sum(w)^2
    vmuhat <- sum(w^2*(y - muhat)^2/(1 - h))/sum(w)^2
  } else if (vartype == "CL") {
    vmuhat <- 1/sum(w)
  } else {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/sum(w)
  }

  lpi <- muhat - stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  res <- list(muhat = muhat, lpi = lpi, upi = upi, 
              lci = lci, uci = uci, tau2h = tau2h, vartype = vartype)

  return(res)

}


#' @export
htsreml <- pima_htsreml
