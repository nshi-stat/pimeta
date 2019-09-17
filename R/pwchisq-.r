#' The Distribution of a Positive Linear Combination of Chiqaure Random Variables
#' 
#' The cumulative distribution function for the distribution of a positive linear
#' combination of \eqn{\chi^2} random variables with the weights
#' (\eqn{\lambda_1, \ldots, \lambda_K}), degrees of freedom
#' (\eqn{\nu_1, \ldots, \nu_K}), and non-centrality parameters
#' (\eqn{\delta_1, \ldots, \delta_K}).
#' 
#' @name pwchisq
#' @rdname pwchisq
#' @param x numeric; value of x > 0 (\eqn{P[X \le x]}).
#' @param lambda numeric vector; weights (\eqn{\lambda_1, \ldots, \lambda_K}).
#' @param nu integer vector; degrees of freedom (\eqn{\nu_1, \ldots, \nu_K}).
#' @param delta numeric vector; non-centrality parameters (\eqn{\delta_1, \ldots, \delta_K}).
#' @param mode numeric; the mode of calculation (see Farabrother, 1984)
#' @param maxit1 integer; the maximum number of iteration.
#' @param eps numeric; the desired level of accuracy.
#' @return
#' \itemize{
#' \item \code{prob}: the distribution function.
#' }
#' @references
#' Farebrother, R. W. (1984).
#' Algorithm AS 204: the distribution of a positive linear combination of
#' \eqn{\chi^2} random variables.
#' \emph{J R Stat Soc Ser C Appl Stat.}
#' \strong{33}(3): 332--339.
<<<<<<< HEAD
#' \url{https://rss.onlinelibrary.wiley.com/doi/10.2307/2347721}.
=======
#' \url{https://www.jstor.org/stable/2347721}.
>>>>>>> ea4600ceb5a8b1e7cc0d7d0d1e90113cff3ab38d
#' @examples
#' # Table 1 of Farabrother (1984)
#' # Q6 (The taget values are 0.0061, 0.5913, and 0.9779)
#' 
#' pimeta::pwchisq( 20, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#' pimeta::pwchisq(100, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#' pimeta::pwchisq(200, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#' # [1] 0.006117973
#' # [1] 0.5913421
#' # [1] 0.9779184
#' @export
pwchisq <- function(x, lambda = 1, nu = 1, delta = 0, mode = 1,
                    maxit1 = 100000, eps = 10^(-10)) {
  
  # initial check
  util_check_nonneg(x)
  util_check_nonneg(lambda)
  util_check_nonneg(nu)
  util_check_nonneg(delta)
  util_check_nonneg(mode)
  util_check_gt(maxit1, 1)
  util_check_nonneg(eps)
  
  if (length(lambda) != length(nu) || length(nu) != length(delta)) {
    stop("'lambda', 'nu', and 'delta' should have the same length.")
  }
  
  # calculate
  res <- pwchisqCpp(
    q      = as.double(x),
    lambda = as.vector(lambda),
    mult   = as.vector(nu),
    delta  = as.vector(delta),
    n      = as.integer(length(lambda)),
    mode   = as.double(mode),
    maxit  = as.integer(maxit1),
    eps    = as.double(eps)
  )
  
  if (res$ifault == -1 || res$ifault == 2) {
    warning("Input value(s) is not appropriate.")
  } else if (res$ifault == 4) {
    warning("The required accuracy could not be achived in 'maxit1' iterations.")
  } else if (res$ifault == 1 || res$ifault == 3 || res$ifault == 5 || res$ifault == 6 ||
             res$ifault == 9 || res$ifault == 10) {
    warning("The estimate of the probability has a problem (e.g., prob < 0, prob > 1, etc.).")
  }
  
  return(res$prob)
  
}