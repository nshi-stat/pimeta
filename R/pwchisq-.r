#' The Distribution of a Positive Linear Combination of Chiqaure Random Variables
#' 
#' The distribution function for the distribution of a positive linear combination
#' of \eqn{\chi^2} random variables with the weights (\eqn{\lambda_1, \ldots, \lambda_K}),
#' degrees of freedom (\eqn{m_1, \ldots, m_K}), and noncentraliy parameters
#' (\eqn{\delta_1, \ldots, \delta_K}).
#' 
#' @name pwchisq
#' @rdname pwchisq
#' @param x numeric; value of x > 0 (\eqn{P[X \leq x]}).
#' @param lambda numeric vector; weights (\eqn{\lambda_1, \ldots, \lambda_K}).
#' @param mult integer vector; degrees of freedom (\eqn{m_1, \ldots, m_K}).
#' @param delta numeric vector; noncentraliy parameters (\eqn{\delta_1, \ldots, \delta_K}).
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
#' \url{https://doi.org/10.2307/2347721}.
#' @examples
#' # Table 1 of Farabrother (1984)
#' # Q1 (The values were 0.0542, 0.4936, 0.8760)
#' pimeta::pwchisq( 1, lambda = c(6,3,1), mult = c(1,1,1), delta = c(0,0,0))
#' pimeta::pwchisq( 7, lambda = c(6,3,1), mult = c(1,1,1), delta = c(0,0,0))
#' pimeta::pwchisq(20, lambda = c(6,3,1), mult = c(1,1,1), delta = c(0,0,0))
#' @export
pwchisq <- function(x, lambda = 1, mult = 1, delta = 0, mode = 1,
                    maxit1 = 100000, eps = 10^(-10)) {
  
  # initial check
  util_check_nonneg(x)
  util_check_nonneg(lambda)
  util_check_nonneg(mult)
  util_check_nonneg(delta)
  util_check_nonneg(mode)
  util_check_nonneg(maxit1)
  util_check_nonneg(eps)
  
  if (length(lambda) != length(mult) || length(mult) != length(delta)) {
    stop("'lambda', 'mult', and 'delta' should have the same length.")
  }
      
  # calculate
  # q, lambda, mult, delta, n, mode, maxit, eps)
  res <- pwchisqCpp(
    q      = as.double(x),
    lambda = as.vector(lambda),
    mult   = as.vector(mult),
    delta  = as.vector(delta),
    n      = as.integer(length(lambda)),
    mode   = as.double(mode),
    maxit  = as.integer(maxit1),
    eps    = as.double(eps)
  )
  
  return(res$prob)
  
}