#' Calculating Heterogeneity Variance
#' 
#' Returns a heterogeneity variance estimate and
#' its confidence interval.
#' 
#' @details Excellent reviews of heterogeneity variance estimation
#' have been published (Sidik & Jonkman, 2007; Veroniki, et al., 2016;
#' Langan, et al., 2018).
#' 
#' @name tau2h
#' @rdname tau2h
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @param maxiter the maximum number of iterations
#' @param method the calculation method for heterogeneity variance (default = "DL").
#' \itemize{
#' \item \code{DL}: DerSimonian--Laird estimator (DerSimonian & Laird, 1986).
#' \item \code{VC}: Variance component type estimator (Hedges, 1983).
#' \item \code{PM}: Paule--Mandel estimator (Paule & Mandel, 1982).
#' \item \code{HM}: Hartung--Makambi estimator (Hartung & Makambi, 2003).
#' \item \code{HS}: Hunter--Schmidt estimator (Hunter & Schmidt, 2004).
#'                  This estimator has negative bias (Viechtbauer, 2005).
#' \item \code{ML}: Maximum likelihood (ML) estimator (e.g., DerSimonian & Laird, 1986).
#' \item \code{REML}: Restricted maximum likelihood (REML) estimator (e.g., DerSimonian & Laird, 1986).
#' \item \code{AREML}: Approximate restricted maximum likelihood estimator (Thompson & Sharp, 1999).
#' \item \code{SJ}: Sidik--Jonkman estimator (Sidik & Jonkman, 2005).
#' \item \code{SJ2}: Sidik--Jonkman improved estimator (Sidik & Jonkman, 2007).
#' \item \code{EB}: Empirical Bayes estimator (Morris, 1983).
#' \item \code{BM}: Bayes modal estimator (Chung, et al., 2013).
#' }
#' @param methodci the calculation method for a confidence interval of
#'                 heterogeneity variance (default = NA).
#' \itemize{
#' \item \code{NA}: a confidence interval will not be calculated.
#' \item \code{ML}: Wald confidence interval with a ML estimator (Biggerstaff & Tweedie, 1997).
#' \item \code{REML}: Wald confidence interval with a REML estimator (Biggerstaff & Tweedie, 1997).
#' }
#' @param alpha the alpha level of the confidence interval
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' \item \code{lci}, \code{uci}: the lower and upper confidence limits
#'       \eqn{\hat{\tau}^2_l} and \eqn{\hat{\tau}^2_u}.
#' }
#' @references
#' Sidik, K., and Jonkman, J. N. (2007).
#' A comparison of heterogeneity variance estimators
#' in combining results of studies.
#' \emph{Stat Med.}
#' \strong{26}(9): 1964-1981.
#' \url{https://doi.org/10.1002/sim.2688}.
#' 
#' Veroniki, A. A., Jackson, D., Viechtbauer, W.,
#' Bender, R., Bowden, J., Knapp, G., Kuss, O.,
#' Higgins, J. P. T., Langan, D., and Salanti, J. (2016).
#' Methods to estimate the between-study variance
#' and its uncertainty in meta-analysis.
#' \emph{Res Syn Meth.}
#' \strong{7}(1): 55-79.
#' \url{https://doi.org/10.1002/jrsm.1164}.
#' 
#' Langan, D., Higgins, J. P. T., Jackson, D.,
#' Bowden, J., Veroniki, A. A., Kontopantelis, E.,
#' Viechtbauer, W., and Simmonds, M. (2018).
#' A comparison of heterogeneity variance estimators in
#' simulated random-effects meta-analyses.
#' \emph{Res Syn Meth.}
#' \emph{In press.}
#' \url{https://doi.org/10.1002/jrsm.1316}.
#' 
#' DerSimonian, R., and Laird, N. (1986).
#' Meta-analysis in clinical trials.
#' \emph{Control Clin Trials.}
#' \strong{7}(3): 177-188.
#' \url{https://doi.org/10.1016/0197-2456(86)90046-2}.
#' 
#' Hedges, L. V. (1983).
#' A random effects model for effect sizes.
#' \emph{Psychol Bull.}
#' \strong{93}(2): 388-395.
#' \url{https://doi.org/10.1037/0033-2909.93.2.388}.
#' 
#' Paule, R. C., and Mandel, K. H. (1982).
#' Consensus values and weighting factors.
#' \emph{J Res Natl Inst Stand Techno.}
#' \strong{87}(5): 377-385.
#' \url{https://doi.org/10.6028/jres.087.022}.
#' 
#' Hartung, J., and Makambi, K. H. (2003).
#' Reducing the number of unjustified significant results
#' in meta-analysis.
#' \emph{Commun Stat Simul Comput.}
#' \strong{32}(4): 1179-1190.
#' \url{https://doi.org/10.1081/SAC-120023884}.
#' 
#' Hunter, J. E., and Schmidt, F. L. (2004).
#' \emph{Methods of Meta-Analysis: Correcting Error and Bias in Research Findings. 2nd edition.}
#' Sage Publications, Inc.
#' 
#' Viechtbauer, W. (2005).
#' Bias and efficiency of meta-analytic variance
#' estimators in the random-effects model.
#' \emph{J Educ Behav Stat.}
#' \strong{30}(3): 261-293.
#' \url{https://doi.org/10.3102/10769986030003261}.
#' 
#' Thompson, S. G., and Sharp, S. J. (1999).
#' Explaining heterogeneity in meta-analysis: a comparison of methods.
#' \emph{Stat Med.}
#' \strong{18}(20): 2693-2708.
#' \url{https://doi.org/10.1002/(SICI)1097-0258(19991030)18:20<2693::AID-SIM235>3.0.CO;2-V}.
#' 
#' Sidik, K., and Jonkman, J. N. (2005).
#' Simple heterogeneity variance estimation for meta-analysis.
#' \emph{J R Stat Soc Ser C Appl Stat.}
#' \strong{54}(2): 367-384.
#' \url{https://doi.org/10.1111/j.1467-9876.2005.00489.x}.
#' 
#' Morris, C. N. (1983).
#' Parametric empirical Bayes inference: theory and applications.
#' \emph{J Am Stat Assoc.}
#' \strong{78}(381): 47-55.
#' \url{https://doi.org/10.1080/01621459.1983.10477920}.
#' 
#' Chung, Y. L., Rabe-Hesketh, S., and Choi, I-H. (2013).
#' Avoiding zero between-study variance estimates
#' in random-effects meta-analysis.
#' \emph{Stat Med.}
#' \strong{32}(23): 4071-4089.
#' \url{https://doi.org/10.1002/sim.5821}.
#' 
#' Biggerstaff, B. J., and Tweedie, R. L. (1997).
#' Incorporating variability in estimates of heterogeneity
#' in the random effects model in meta-analysis.
#' \emph{Stat Med.}
#' \strong{16}(7): 753-768.
#' \url{https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<753::AID-SIM494>3.0.CO;2-G}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h(sbp$y, sbp$sigmak)
#' @export
tau2h <- function(y, se, maxiter = 100, method = c("DL", "VC", "PM", "HM", "HS", "ML", "REML",
                  "AREML", "SJ", "SJ2", "EB", "BM"), methodci = c(NA, "ML", "REML"), alpha = 0.05) {
  
  # initial check
  lstm <- c("DL", "VC", "PM", "HM", "HS", "ML", "REML", "AREML", "SJ", "SJ2", "EB", "BM")
  lstc <- c(NA, "ML", "REML")
  method <- match.arg(method)
  methodci <- match.arg(methodci)
  
  util_check_num(y)
  util_check_nonneg(se)
  util_check_gt(maxiter, 1)
  util_check_inrange(alpha, 0.0, 1.0)

  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (!is.element(method, lstm)) {
    stop("Unknown 'method' specified.")
  } else if (!is.element(methodci, lstc)) {
    stop("Unknown 'methodci' specified.")
  }

  # point estiamte
  if (method == "DL") {
    res <- tau2h_dl(y = y, se = se)
  } else if (method == "VC") {
    res <- tau2h_vc(y = y, se = se)
  } else if (method == "PM") {
    res <- tau2h_pm(y = y, se = se)
  } else if (method == "HM") {
    res <- tau2h_hm(y = y, se = se)
  } else if (method == "HS") {
    res <- tau2h_hs(y = y, se = se)
  } else if (method == "ML") {
    res <- tau2h_ml(y = y, se = se, maxiter = maxiter)
  } else if (method == "REML") {
    res <- tau2h_reml(y = y, se = se, maxiter = maxiter)
  } else if (method == "AREML") {
    res <- tau2h_areml(y = y, se = se, maxiter = maxiter)
  } else if (method == "SJ") {
    res <- tau2h_sj(y = y, se = se)
  } else if (method == "SJ2") {
    res <- tau2h_sj2(y = y, se = se)
  } else if (method == "EB") {
    res <- tau2h_eb(y = y, se = se, maxiter = maxiter)
  } else if (method == "BM") {
    res <- tau2h_bm(y = y, se = se, maxiter = maxiter)
  }

  # confidence interval
  if (is.na(methodci)) {
  } else if (methodci == "ML") {
    tau2hml <- tau2h_ml(y = y, se = se, maxiter = maxiter)
    resci <- tau2h_wald_ml(se = se, tau2h = tau2hml$tau2h, alpha = alpha)
    res <- append(res, resci)
  } else if (methodci == "REML") {
    tau2hreml <- tau2h_reml(y = y, se = se, maxiter = maxiter)
    resci <- tau2h_wald_reml(se = se, tau2h = tau2hreml$tau2h, alpha = alpha)
    res <- append(res, resci)
  }
  
  res <- append(res, list(y = y, se = se, method = method, methodci = methodci))
  class(res) <- "pima_tau2h"
  
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
#' @method print pima_tau2h
print.pima_tau2h <- function(x, digits = 3, ...) {
  
  cat("\nHeterogeneity Variance for Random-Effects Meta-Analysis\n\n")
  
  cat("Estimation method(s)\n")
  if (x$method == "DL") {
    cat(" Point estimation: Dersimonian-Laird estimator\n")
  } else if (x$method == "VC") {
    cat(" Point estimation: Variance component type estimator\n")
  } else if (x$method == "PM") {
    cat(" Point estimation: Paule--Mandel estimator\n")
  } else if (x$method == "HM") {
    cat(" Point estimation: Hartung--Makambi estimator\n")
  } else if (x$method == "HS") {
    cat(" Point estimation: Hunter--Schmidt estimator\n")
    cat(" # Note: this estimator has negative bias.\n")
  } else if (x$method == "ML") {
    cat(" Point estimation: Maximum likelihood estimator\n")
  } else if (x$method == "REML") {
    cat(" Point estimation: Restricted maximum likelihood estimator\n")
  } else if (x$method == "AREML") {
    cat(" Point estimation: Approximate restricted maximum likelihood estimator\n")
  } else if (x$method == "SJ") {
    cat(" Point estimation: Sidik--Jonkman estimator\n")
  } else if (x$method == "SJ2") {
    cat(" Point estimation: Sidik--Jonkman improved estimator\n")
  } else if (x$method == "EB") {
    cat(" Point estimation: Empirical Bayes estimator\n")
  } else if (x$method == "BM") {
    cat(" Point estimation: Bayes modal estimator\n")
  }
  if (is.na(x$methodci)) {
  } else if (x$methodci == "ML") {
    cat(" Confidence interval: Maximum likelihood estimator\n")
  } else if (x$methodci == "REML") {
    cat(" Confidence interval: Restricted maximum likelihood estimator\n")
  }
  cat("\n")
  
  cat(paste0("No. of studies: ", length(x$y), "\n\n"))
  
  if (is.na(x$methodci)) {
    cat(paste0("tau-squared: ", format(round(x$tau2h, digits), nsmall = digits), "\n"))
  } else {
    cat(paste0("tau-squared [", (1 - x$alpha)*100, "% confidence interval]: "))
    cat(paste0(format(round(x$tau2h, digits), nsmall = digits), " [",
               format(round(x$lci, digits), nsmall = digits), ", ",
               format(round(x$uci, digits), nsmall = digits), "]\n"))
  }
  cat("\n")

  invisible(x)
  
}
