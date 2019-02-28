#' Partlett--Riley prediction interval
#'
#' A function for the Partlett--Riley PIs
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
#' \item \code{HK}: the Hartung (1999)'s estimator (the Hartung and Knapp (2001)'s estimator).
#' \item \code{SJBC}: the Sidik and Jonkman (2006)'s bias coreccted estimator.
#' \item \code{KR}: the Kenward and Roger (1997)'s approach.
#' \item \code{APX}: an approximate estimator (DerSimonian & Laird, 1986), \eqn{(1/\sum{\hat{w}_i})^{-1}}.
#' }
#' @param maxiter the maximum number of iterations
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' \item \code{vmuhat}: the variance estimate for \eqn{\hat{\mu}}.
#' \item \code{nup}: degrees of freedom for the prediction interval.
#' \item \code{nuc}: degrees of freedom for the confidence interval.
#' }
#' @references
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
#' \url{https://doi.org/10.1002/sim.7140}
#' 
#' Hartung, J. (1999).
#' An alternative method for meta-analysis.
#' \emph{Biom J.}
#' \strong{41}(8): 901-916.
#' \url{https://doi.org/10.1002/(SICI)1521-4036(199912)41:8<901::AID-BIMJ901>3.0.CO;2-W}
#' 
#' Hartung, J., and Knapp, G. (2001).
#' On tests of the overall treatment effect in meta-analysis with
#' normally distributed responses.
#' \emph{Stat Med.}
#' \strong{20}(12): 1771-1782.
#' \url{https://doi.org/10.1002/sim.791}
#' 
#' Sidik, K., and Jonkman, J. N. (2006).
#' Robust variance estimation for random effects meta-analysis.
#' \emph{Comput Stat Data Anal.}
#' \strong{50}(12): 3681-3701.
#' \url{https://doi.org/10.1016/j.csda.2005.07.019}
#' 
#' Kenward, M. G., and Roger, J. H. (1997).
#' Small sample inference for fixed effects from restricted
#' maximum likelihood.
#' \emph{Biometrics.}
#' \strong{53}(3): 983-997.
#' \url{https://doi.org/10.2307/2533558}
#' 
#' DerSimonian, R., and Laird, N. (1986).
#' Meta-analysis in clinical trials.
#' \emph{Control Clin Trials.}
#' \strong{7}(3): 177-188.
#' @seealso
#' \code{\link[=pima]{pima}}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::pima_htsreml(sbp$y, sbp$sigmak)
#' @export
pima_htsreml <- function(y, sigma, alpha = 0.05,
                         vartype = c("HK", "SJBC", "KR", "CL", "APX"), maxiter = 100) {

  # initial check
  lstm <- c("HK", "SJBC", "KR", "CL", "APX")
  vartype <- match.arg(vartype)
  
  util_check_num(y)
  util_check_num(sigma)
  util_check_num(alpha)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (!is.element(vartype, lstm)) {
    stop("Unknown 'vartype' specified.")
  }
  
  # estimation
  k <- length(y)
  tau2h <- tau2h_reml(y = y, se = sigma, maxiter = maxiter)$tau2h
  
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  w1p <- sum(w)
  w2p <- sum(w^2)
  w3p <- sum(w^3)
  nup <- k - 2
  nuc <- k - 1
  
  if (vartype == "HK") {
    q <- sum(w*(y - muhat)^2)/(k - 1.0)
    vmuhat <- max(1.0, q)/w1p
    method <- "HK"
  } else if (vartype == "SJBC") {
    lambda <- sigma^2 + tau2h
    h <- 2.0*w/w1p - sum(w^2*lambda)/(lambda*w1p^2)
    vmuhat <- sum(w^2*(y - muhat)^2/(1.0 - h))/w1p^2
    method <- "SJ"
  } else if (vartype == "KR") {
    IE <- w2p*0.5 - w3p/w1p + 0.5*(w2p/w1p)^2
    vmuhat <- (1.0 + 2.0*(w3p/w1p - (w2p/w1p)^2)/IE)/w1p
    nup <- 2*IE/(vmuhat*w2p)^2 - 1
    nuc <- 2*IE/(vmuhat*w2p)^2
    method <- "KR"
  } else if (vartype == "CL" | vartype == "APX") {
    vmuhat <- 1.0/w1p
    method <- "APX"
  } else {
    vmuhat <- sum(w*(y - muhat)^2)/(k - 1)/w1p
    method <- "HK"
  }

  lpi <- muhat - stats::qt(1.0 - alpha*0.5, nup)*sqrt(tau2h + vmuhat)
  upi <- muhat + stats::qt(1.0 - alpha*0.5, nup)*sqrt(tau2h + vmuhat)
  lci <- muhat - stats::qt(1.0 - alpha*0.5, nuc)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1.0 - alpha*0.5, nuc)*sqrt(vmuhat)
  
  res <- list(muhat = muhat, lpi = lpi, upi = upi, lci = lci, uci = uci, 
              tau2h = tau2h, vmuhat = vmuhat, nup = nup, nuc = nuc,
              method = method, y = y, se = sigma, alpha = alpha)

  return(res)

}


#' @export
htsreml <- pima_htsreml
