#' A parametric bootstrap prediction interval
#' 
#' A subroutine for the parametric bootstrap PI
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
#' @seealso
#' \code{\link[=pima]{pima()}}.
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::pima_boot(sbp$y, sbp$sigmak, B = 50000)}
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # A parametric bootstrap prediction interval
#' #  Heterogeneity variance: DerSimonian-Laird
#' #  SE for average treatment effect: Hartung
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3341 [-0.8769, 0.2248]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3341 [-0.5660, -0.0976]
#' # 
#' # Heterogeneity variance (tau^2):
#' #  0.0282
#' # 
#' @export
pima_boot <- function(y, sigma, alpha = 0.05, B = 25000, maxit1 = 100000,
                   eps = 10^(-10), lower = 0, upper = 1000, maxit2 = 1000,
                   tol = .Machine$double.eps^0.25, rnd = NULL) {

  ## .. need more more strictry check.
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  } else if (min(sigma) < 0.0) {
    stop("'sigma' should be positive.")
  } else if (B < 1) {
    stop("'B' should be grater than 1.")
  }

  if (is.null(rnd)) {
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
  } else {
    rndtau2 <- rnd
  }

  tau2h <- tau2h_dl(y = y, se = sigma)
  w <- (sigma^2 + tau2h)^-1
  muhat <- list(muhat = sum(y*w) / sum(w))
  res <- bootPICppWrap(
    rnd   = as.vector(rndtau2),
    y     = as.vector(y),
    sigma = as.vector(sigma),
    alpha = as.double(alpha)
  )
  
  res <- append(append(muhat, res),
                list(tau2h = tau2h, method = "boot", y = y, se = sigma,
                     alpha = alpha, rnd = rndtau2))

  return(res)

}


#' @export
bootPI <- pima_boot
