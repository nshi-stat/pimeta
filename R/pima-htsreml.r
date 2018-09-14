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
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
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
    method <- "HK"
  } else if (vartype == "SJBC") {
    lambda <- sigma^2 + tau2h
    h <- 2*w/sum(w) - sum(w^2*lambda)/lambda/sum(w)^2
    vmuhat <- sum(w^2*(y - muhat)^2/(1 - h))/sum(w)^2
    method <- "SJ"
  } else if (vartype == "CL") {
    vmuhat <- 1/sum(w)
    method <- "CL"
  } else {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/sum(w)
    method <- "HK"
  }

  lpi <- muhat - stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1-0.05/2, k - 2)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1-0.05/2, k - 1)*sqrt(vmuhat)
  
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci,
              tau2h = tau2h, method = method, y = y, se = sigma, alpha = alpha)
  class(res) <- "pima" 
  
  return(res)

}


#' @export
htsreml <- pima_htsreml
