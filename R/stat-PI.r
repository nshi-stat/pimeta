#' A parametric bootstrap prediction interval (Nagashima et al., 2018)
#'
#' @name bootPI
#' @rdname bootPI
#' @aliases bootPI
#' @param y the effect size estimates vector
#' @param sigma the within studies variances vector
#' @param alpha the alpha level of the prediction interval
#' @param B the number of bootstrap samples
#' @param maxit1 the maximum number of iteration for the exact distribution function of \eqn{Q}
#' @param eps the desired level of accuracy for the exact distribution function of \eqn{Q}
#' @param lower the lower limit of random numbers of \eqn{\tau^2}
#' @param upper the lower upper of random numbers of \eqn{\tau^2}
#' @param maxit2 the maximum number of iteration for numerical inversions
#' @param tol the desired level of accuracy for numerical inversions
#' @return The average treatment effect estimate \eqn{\hat{\mu}} (\code{muhat}),
#' and the lower and upper prediction limits \eqn{\hat{c}_l} (\code{lbpi}) and \eqn{\hat{c}_u} (\code{ubpi}).
#' @references
#' Nagashima, K., Noma, H., and Furukawa, T. A. (2018).
#' Prediction intervals for random-effects meta-analysis:
#' a confidence distribution approach.
#' \emph{Stat Methods Med Res}.
#' \emph{In press}.
#' \url{https://arxiv.org/abs/1804.01054}.
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::bootPI(sbp$y, sbp$sigmak, B = 50000)}
#' # $muhat
#' # [1] -0.3340597
#' # $lbpi
#' # [1] -0.8768976
#' # $ubpi
#' # [1] 0.2248231
#' @export
bootPI <- function(y, sigma, alpha = 0.05, B = 25000, maxit1 = 100000,
                   eps = 10^(-10), lower = 0, upper = 1000, maxit2 = 1000,
                   tol = .Machine$double.eps^0.25) {

  ## .. need more more strictry check.
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (min(sigma) < 0.0) {
    stop("'sigma' should be positive.")
  } else if (B < 1) {
    stop("'B' should be grater than 1.")
  }

  rndtau2 <- rtau2CppWrap(
    n      = as.integer(B),
    y      = as.vector(y),
    sigma  = as.vector(sigma),
    mode   = as.double(1),
    maxit1 = as.integer(maxit1),
    eps    = as.double(eps),
    lower  = as.double(lower),
    upper  = as.double(upper),
    maxit2 = as.integer(maxit2),
    tol    = as.double(tol)
  )

  k <- length(y)
  tauh2 <- max(0, (sum(sigma^-2 * (y - sum(sigma^-2*y) / sum(sigma^-2))^2) - (k - 1)) /
                 (sum(sigma^-2) - sum(sigma^-4)/sum(sigma^-2)))
  w <- (sigma^2 + tauh2)^-1
  muhat <- list(muhat = sum(y*w) / sum(w))
  res <- bootPICppWrap(
    rnd   = as.vector(rndtau2),
    y     = as.vector(y),
    sigma = as.vector(sigma),
    alpha = as.double(alpha)
  )

  return(append(muhat, res))

}


#' Higgins-Thompson-Spiegelhalter prediction interval with the Dersimonian-Laird estimator for \eqn{\hat{\tau}}
#'
#' @name htsdl
#' @rdname htsdl
#' @aliases htsdl
#' @param y the effect size estimates vector
#' @param sigma the within studies variances vector
#' @param alpha the alpha level of the prediction interval
#' @return The average treatment effect estimate \eqn{\hat{\mu}} (\code{muhat}),
#' the lower and upper prediction limits \eqn{\hat{c}_l} (\code{lpi}) and \eqn{\hat{c}_u} (\code{upi}),
#' the DL estimator for \eqn{\hat{\tau}} (\code{tau2h}).
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::htsdl(sbp$y, sbp$sigmak)
#' # $muhat
#' # [1] -0.3340597
#' # $lbpi
#' # [1] -0.7597777
#' # $ubpi
#' # [1] 0.09165839
#' # $tau2h
#' # [1] 0.02824971
#' @export
htsdl <- function(y, sigma, alpha = 0.05) {

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
  res <- list(muhat = muhat, lpi = lpi, upi = upi, tau2h = tau2h)

  return(res)

}


#' Higgins-Thompson-Spiegelhalter prediction interval with a REML estimator for \eqn{\hat{\tau}}
#'
#' @name htsreml
#' @rdname htsreml
#' @aliases htsreml
#' @param y the effect size estimates vector
#' @param sigma the within studies variances vector
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
htsreml <- function(y, sigma, alpha = 0.05, vartype = c("HK", "SJBC", "CL"), maxiter = 100) {

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
  res <- list(muhat = muhat, lpi = lpi, upi = upi, tau2h = tau2h, vartype = vartype)

  return(res)

}
