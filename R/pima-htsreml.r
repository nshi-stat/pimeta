#' Partlett--Riley prediction interval
#'
#' A subroutine for the Partlett--Riley PI
#' based on the REML estimator (Partlett & Riley, 2017)
#' 
#' @name pima_htsreml
#' @rdname pima_htsreml
#' @aliases htsreml
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @param vartype the type of the variance estimator for \eqn{\hat{\mu}} (default = "HK"):
#' \itemize{
#' \item \code{HK}: the Hartung and Knapp (2001)'s estimator.
#' \item \code{SJBC}: the Sidik and Jonkman (2006)'s bias coreccted estimator.
#' \item \code{CL}: a standard estimator, \eqn{(1/\sum{\hat{w}_i})^{-1}}.
#' }
#' @param maxiter the maximum number of iterations
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
#' 
#' Hartung, J., and Knapp, G. (2001).
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Stat Med.}
#' \strong{20}(12): 1771-1782.
#' 
#' Sidik, K., and Jonkman, J. N. (2006).
#' Robust variance estimation for random effects meta-analysis.
#' \emph{Comput Stat Data Anal.}
#' \strong{50}(12): 3681-3701.
#' @seealso
#' \code{\link[=pima]{pima()}}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::pima_htsreml(sbp$y, sbp$sigmak)
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # Partlett-Riley prediction interval
#' #  Heterogeneity variance: REML
#' #  SE for average treatment effect: Hartung-Knapp
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3287 [-0.9887, 0.3312]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3287 [-0.5761, -0.0814]
#' # 
#' # Heterogeneity variance (tau^2):
#' #  0.0700
#' #
#' @export
pima_htsreml <- function(y, sigma, alpha = 0.05,
                         vartype = c("HK", "SJBC", "CL"), maxiter = 100) {

  # initial check
  lstm <- c("HK", "SJBC", "CL")
  method <- match.arg(method)
  
  util_check_num(y)
  util_check_num(sigma)
  util_check_num(alpha)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (!is.element(method, lstm)) {
    stop("Unknown 'method' specified.")
  }
  
  # estimation
  k <- length(y)
  tau2h <- tau2h_reml(y = y, se = sigma, maxiter = maxiter)$tau2h
  
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)

  if (vartype == "HK") {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/sum(w)
    method <- "HK"
  } else if (vartype == "SJBC") {
    lambda <- sigma^2 + tau2h
    h <- 2.0*w/sum(w) - sum(w^2*lambda)/lambda/sum(w)^2
    vmuhat <- sum(w^2*(y - muhat)^2/(1.0 - h))/sum(w)^2
    method <- "SJ"
  } else if (vartype == "CL") {
    vmuhat <- 1/sum(w)
    method <- "CL"
  } else {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/sum(w)
    method <- "HK"
  }

  lpi <- muhat - stats::qt(1.0 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1.0 - alpha*0.5, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1.0 - alpha*0.5, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1.0 - alpha*0.5, k - 1)*sqrt(vmuhat)
  
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci,
              tau2h = tau2h, method = method, y = y, se = sigma, alpha = alpha)

  return(res)

}


#' @export
htsreml <- pima_htsreml
