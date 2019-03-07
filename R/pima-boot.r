#' A parametric bootstrap prediction interval
#' 
#' A function for the parametric bootstrap PI
#' based on confidence distribution (Nagashima et al., 2018).
#' A parametric bootstrap confidence interval is also calculated
#' based on the same sampling method for bootstrap PI.
#'
#' @name pima_boot
#' @rdname pima_boot
#' @aliases bootPI
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @param B the number of bootstrap samples
#' @param maxit1 the maximum number of iteration for the exact distribution function of \eqn{Q}
#' @param eps the desired level of accuracy for the exact distribution function of \eqn{Q}
#' @param lower the lower limit of random numbers of \eqn{\tau^2}
#' @param upper the upper limit of random numbers of \eqn{\tau^2}
#' @param maxit2 the maximum number of iteration for numerical inversions
#' @param tol the desired level of accuracy for numerical inversions
#' @param rnd a vector of random numbers from the exact distribution of \eqn{\tau^2}
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{lpi}, \code{lpi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Nagashima, K., Noma, H., and Furukawa, T. A. (2018).
#' Prediction intervals for random-effects meta-analysis:
#' a confidence distribution approach.
#' \emph{Stat Methods Med Res}.
#' \emph{In press}.
#' \url{https://doi.org/10.1177/0962280218773520}.
#' 
#' Hartung, J. (1999).
#' An alternative method for meta-analysis.
#' \emph{Biom J.}
#' \strong{41}(8): 901-916.
#' \url{https://doi.org/10.1002/(SICI)1521-4036(199912)41:8<901::AID-BIMJ901>3.0.CO;2-W}
#' @seealso
#' \code{\link[=pima]{pima}}.
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::pima_boot(sbp$y, sbp$sigmak, B = 50000)}
#' @export
pima_boot <- function(y, sigma, alpha = 0.05, B = 25000, maxit1 = 100000,
                      eps = 10^(-10), lower = 0, upper = 1000, maxit2 = 1000,
                      tol = .Machine$double.eps^0.25, rnd = NULL, parallel = FALSE,
                      seed = NULL) {
  
  # initial check
  util_check_num(y)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  util_check_gt(B, 1)
  util_check_nonneg(parallel)
  util_check_num(seed)
  util_check_gt(maxit1, 1)
  util_check_gt(eps, 0)
  util_check_ge(lower, 0)
  util_check_gt(upper, 0)
  util_check_gt(maxit2, 1)
  util_check_gt(tol, 0)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (lower >= upper) {
    stop("'upper' should be greater than 'lower'.")
  }
  
  if (B < 5000) {
    warning("'B' > 5000 is recommended.")
  }
  
  k <- length(y)
  
  # random numbers generation
  if (is.null(rnd)) {
    if (parallel == FALSE) {
<<<<<<< HEAD
      parallel <- 1
=======
      parallel = 1
>>>>>>> 0a112d94f05a3b66b0a66c33d0e72212def7faf4
    }
    set.seed(seed)
    rndtau2 <- rtau2CppWrap(
      n       = as.integer(B),
      y       = as.vector(y),
      sigma   = as.vector(sigma),
      mode    = as.double(1),
      maxit1  = as.integer(maxit1),
      eps     = as.double(eps),
      lower   = as.double(lower),
      upper   = as.double(upper),
      maxit2  = as.integer(maxit2),
      tol     = as.double(tol),
      nthread = as.integer(parallel)
    )
  } else {
    rndtau2 <- rnd
  }
  
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- list(muhat = sum(y*w) / sum(w))
  res <- bootPICppWrap(
    rnd   = as.vector(rndtau2),
    y     = as.vector(y),
    sigma = as.vector(sigma),
    alpha = as.double(alpha)
  )
  
  res <- append(append(muhat, res),
                list(tau2h = tau2h, vmuhat = NULL, nup = k - 1, nuc = k - 1,
                     method = "boot", y = y, se = sigma,
                     alpha = alpha, rnd = rndtau2))
  
  return(res)
  
}


#' @export
bootPI <- pima_boot
