#' Calculating Prediction Intervals
#' 
#' This function calculates prediction intervals using \code{pima_boot},
#' \code{pima_hts}, or \code{pima_htsreml} functions.
#'
#' @name pima
#' @rdname pima
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @param method the calculation method for the pretiction interval (default = "boot").
#' \itemize{
#' \item \code{boot}: A parametric bootstrap prediction interval.
#' \item \code{HTS}: the Higgins--Thompson--Spiegelhalter prediction interval / 
#'            the DerSimonian & Laird estimator for \eqn{\tau^2} with
#'            a standard SE estimator for the average effect,
#'            \eqn{(1/\sum{\hat{w}_i})^{-1}}.
#' \item \code{HK}: Partlett & Riley's prediction interval /
#'            the REML estimator for \eqn{\tau^2} with
#'            the Hartung and Knapp (2001)'s SE estimator for the average effect.
#' \item \code{SJ}: Partlett & Riley's prediction interval /
#'            the REML estimator for \eqn{\tau^2} with
#'            the Sidik and Jonkman (2006)'s bias coreccted SE estimator
#'             for the average effect.
#' \item \code{CL}: a prediction interval with REML and standard SE /
#'            the REML estimator for \eqn{\tau^2} with
#'            a standard SE estimator for the average effect.
#' }
#' @param B the number of bootstrap samples
#' @param maxit1 the maximum number of iteration for the exact distribution function of \eqn{Q}
#' @param eps the desired level of accuracy for the exact distribution function of \eqn{Q}
#' @param lower the lower limit of random numbers of \eqn{\tau^2}
#' @param upper the lower upper of random numbers of \eqn{\tau^2}
#' @param maxit2 the maximum number of iteration for numerical inversions
#' @param tol the desired level of accuracy for numerical inversions
#' @param rnd a vector of random numbers from the exact distribution of \eqn{\tau^2}
#' @param maxiter the maximum number of iteration for REML estimation
#' @return The average treatment effect estimate \eqn{\hat{\mu}} (\code{muhat}),
#' and the lower and upper prediction limits \eqn{\hat{c}_l} (\code{lbpi}) and \eqn{\hat{c}_u} (\code{ubpi}).
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' 
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
#' 
#' Nagashima, K., Noma, H., and Furukawa, T. A. (2018).
#' Prediction intervals for random-effects meta-analysis:
#' a confidence distribution approach.
#' \emph{Stat Methods Med Res}.
#' \emph{In press}.
#' \url{https://doi.org/10.1177/0962280218773520}.
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, B = 50000)}
#' # Average treatment effect [95%PI]:
#' # -0.3341 [-0.8769, 0.2248]
#' # Average treatment effect [95%CI]:
#' # -0.3341 [-0.5660, -0.0976]
#' # Heterogeneity variance (tau^2):
#' # 0.0282
#'
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "HTS")}
#' # Average treatment effect [95%PI]:
#' # -0.3341 [-0.7598, 0.0917]
#' # Average treatment effect [95%CI]:
#' # -0.3341 [-0.5068, -0.1613]
#' # Heterogeneity variance (tau^2):
#' # 0.0282
#' 
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "HK")}
#' # Average treatment effect [95%PI]:
#' # -0.3287 [-0.9887, 0.3312]
#' # Average treatment effect [95%CI]:
#' # -0.3287 [-0.5761, -0.0814]
#' # Heterogeneity variance (tau^2):
#' # 0.0700
#' 
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "SJ")}
#' # Average treatment effect [95%PI]:
#' # -0.3287 [-0.9835, 0.3261]
#' # Average treatment effect [95%CI]:
#' # -0.3287 [-0.5625, -0.0950]
#' # Heterogeneity variance (tau^2):
#' # 0.0700
#' @export
pima <- function(y, se, alpha = 0.05, method = "boot",
                 B = 25000, maxit1 = 100000, eps = 10^(-10), lower = 0, upper = 1000,
                 maxit2 = 1000, tol = .Machine$double.eps^0.25, rnd = NULL,
                 maxiter = 100) {
  
  ## .. may be need more more strictry check.
  
  if (!is.element(method, c("boot", "HTS", "HK", "SJ", "CL")))
    stop("Unknown 'method' specified.")
  
  if (length(se) != length(y)) 
    stop("'y' and 'se' should have the same length.")
  
  if (min(se) < 0.0)
    stop("'se' should be positive.")
  
  if (B < 1)
    stop("'B' should be grater than 1.")
  
  
  if (method == "boot") {
    res <- pima_boot(y      = y, 
                     sigma  = se, 
                     alpha  = alpha,
                     B      = B,
                     maxit1 = maxit1,
                     eps    = eps, 
                     lower  = lower,
                     upper  = upper, 
                     maxit2 = maxit2,
                     rnd    = rnd)
  } else if (method == "HTS") {
    res <- pima_hts(y      = y, 
                    sigma  = se, 
                    alpha  = alpha)
  } else if (method == "HK") {
    res <- pima_htsreml(y       = y, 
                        sigma   = se, 
                        alpha   = alpha,
                        vartype = "HK",
                        maxiter = maxiter)
  } else if (method == "SJ") {
    res <- pima_htsreml(y       = y, 
                        sigma   = se, 
                        alpha   = alpha,
                        vartype = "SJBC",
                        maxiter = maxiter)
  } else if (method == "CL") {
    res <- pima_htsreml(y       = y, 
                        sigma   = se, 
                        alpha   = alpha,
                        vartype = "CL",
                        maxiter = maxiter)
  }
  
  res <- append(res, list(method = method, y = y, se = se, alpha = alpha))
  class(res) <- "pima" 
  
  return(res)
  
}


#' Print Results
#' 
#' \code{print} prints its argument and returns it invisibly (via \code{invisible(x)}).
#'
#' @param x print to display
#' @param digits a value for digits specifies the minimum number
#'               of significant digits to be printed in values.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method print pima
print.pima <- function(x, digits = 4, ...) {
  
  cat("\nPrediction Interval for Random-Effects Meta-Analysis\n\n")
  
  if (x$method == "boot") {
    cat(paste0("A parametric bootstrap prediction interval\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: Hartung\n\n"))
  } else if (x$method == "HTS") {
    cat(paste0("Higgins-Thompson-Spiegelhalter prediction interval\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: standard\n\n"))
  } else if (x$method == "HK") {
    cat(paste0("Partlett-Riley prediction interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Hartung-Knapp\n\n"))
  } else if (x$method == "SJ") {
    cat(paste0("Partlett-Riley prediction interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Sidik-Jonkman\n\n"))
  } else if (x$method == "CL") {
    cat(paste0("A prediction interval with REML and standard SE\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: standard\n\n"))
  }  
  cat(paste0("Average treatment effect [", (1 - x$alpha)*100, "%PI]:\n"))
  cat(paste0(" ", format(round(x$muhat, digits), nsmall = digits), " [",
             format(round(x$lpi, digits), nsmall = digits), ", ",
             format(round(x$upi, digits), nsmall = digits), "]\n\n"))
  
  cat(paste0("Average treatment effect [", (1 - x$alpha)*100, "%CI]:\n"))
  cat(paste0(" ", format(round(x$muhat, digits), nsmall = digits), " [",
             format(round(x$lci, digits), nsmall = digits), ", ",
             format(round(x$uci, digits), nsmall = digits), "]\n\n"))
  
  cat(paste0("Heterogeneity variance (tau^2):\n"))
  cat(paste0(" ", format(round(x$tau2, digits), nsmall = digits), "\n\n"))
  
  invisible(x)
  
}