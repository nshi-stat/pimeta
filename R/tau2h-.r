#' Calculating Heterogeneity Variance
#' 
#' Returns a heterogeneity variance estimate and
#' its confidence interval.
#' 
#' @details Excellent reviews of heterogeneity variance estimation
#' have been published (Sidik & Jonkman, 2007; Veroniki, et al., 2016).
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
#' \item \code{HM}: Hartung--Makambi estimator (Hartung & Makambi, 1983).
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
#'                 heterogeneity variance (default = "ML").
#' \itemize{
#' \item \code{ML}: Wald confidence interval with a ML estimator (Biggerstaff & Tweedie, 1997).
#' }
#' @param alpha the alpha level of the confidence interval
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' \item \code{lci}, \code{uci}: the lower and upper confidence limits
#'       \eqn{\hat{\tau}^2_l} and \eqn{\hat{\tau}^2_u}.
#' }
#' @references
#' Veroniki, A. A., Jackson, D., Viechtbauer, W.,
#' Bender, R., Bowden, J., Knapp, G., Kuss, O.,
#' Higgins, J. P. T., Langan, D., and Salanti, J. (2016).
#' Methods to estimate the between-study variance
#' and its uncertainty in meta-analysis.
#' \emph{Res Syn Meth.}
#' \strong{7}(1): 55-79.
#' 
#' Sidik, K., and Jonkman, J. N. (2007).
#' A comparison of heterogeneity variance estimators
#' in combining results of studies.
#' \emph{Stat Med.}
#' \strong{26}(9): 1964-1981.
#' 
#' DerSimonian, R., and Laird, N. (1986).
#' Meta-analysis in clinical trials.
#' \emph{Control Clin Trials.}
#' \strong{7}(3): 177-188.
#' 
#' Hedges, L. V. (1983).
#' A random effects model for effect sizes.
#' \emph{Psychol Bull.}
#' \strong{93}(2): 388-395.
#' 
#' Paule, R. C., and Mandel, K. H. (1982).
#' Consensus values and weighting factors.
#' \emph{J Res Natl Inst Stand Techno.}
#' \strong{87}(5): 377-385.
#' 
#' Hartung, J., and Makambi, K. H. (2003).
#' Reducing the number of unjustified significant results
#' in meta-analysis.
#' \emph{Commun Stat Simul Comput.}
#' \strong{32}(4): 1179-1190.
#' 
#' Hunter, J. E., and Schmidt, F. L. (2004).
#' \emph{Methods of Meta-Analysis: Correcting Error and Bias in Research Findings. 2nd edition.}
#' Sage Publications, Inc.
#' 
#' Viechtbauer, W. (2005).
#' Bias and efficiency of meta-analytic variance
#' estimators in the random-effects Model.
#' \emph{J Educ Behav Stat.}
#' \strong{30}(3): 261-293.
#' 
#' Thompson, S. G., and Sharp, S. J. (1999).
#' Explaining heterogeneity in meta-analysis: a comparison of methods.
#' \emph{Stat Med.}
#' \strong{18}(20): 2693-2708.
#' 
#' Sidik, K., and Jonkman, J. N. (2005).
#' Simple heterogeneity variance estimation for meta-analysis.
#' \emph{J R Stat Soc Ser C Appl Stat.}
#' \strong{54}(2): 367-384.
#' 
#' Morris, C. N. (1983).
#' Parametric empirical Bayes inference: theory and applications.
#' \emph{J Am Stat Assoc.}
#' \strong{78}(381): 47-55.
#' 
#' Chung, Y. L., Rabe-Hesketh, S., and Choi, I-H. (2013).
#' Avoiding zero between-study variance estimates
#' in random-effects meta-analysis.
#' \emph{Stat Med.}
#' \strong{32}(23): 4071-4089.
#' 
#' Biggerstaff, B. J., and Tweedie, R. L. (1997).
#' Incorporating variability in estimates of heterogeneity
#' in the random effects model in meta-analysis.
#' \emph{Stat Med.}
#' \strong{16}(7): 753-768.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h(sbp$y, sbp$sigmak, method = "DL", methodci = "ML")
#' @export
tau2h <- function(y, se, maxiter = 100, method = c("DL", "VC", "PM", "HM", "HS", "ML", "REML",
                  "AREML", "SJ", "SJ2", "EB", "BM"), methodci = c("ML", "REML"), alpha = 0.05) {
  
  # initial check
  lstm <- c("DL", "VC", "PM", "HM", "HS", "ML", "REML", "AREML", "SJ", "SJ2", "EB", "BM")
  lstc <- c("ML", "REML")
  method <- match.arg(method)
  methodci <- match.arg(methodci)
  
  if (is.null(y)) {
    stop("'y' is a null value.")
  } else if (is.null(se)) {
    stop("'se' is a null value.")
  } else if (any(is.na(y))) {
    stop("'y' has missing value(s).")
  } else if (any(is.na(se))) {
    stop("'se' has missing value(s).")
  } else if (any(is.infinite(y))) {
    stop("'y' has infinite value(s).")
  } else if (any(is.infinite(se))) {
    stop("'se' has infinite value(s).")
  } else if (any(is.nan(y))) {
    stop("'y' has NaN(s).")
  } else if (any(is.nan(se))) {
    stop("'se' has NaN(s).")
  } else if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  } else if (!is.element(method, lstm)) {
    stop("Unknown 'method' specified.")
  } else if (!is.element(methodci, lstc)) {
    stop("Unknown 'methodci' specified.")
  } else if (alpha < 0.0 | alpha > 1.0) {
    stop("'alpha' should be 0.0 < alpha < 1.0.")
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
  if (methodci == "ML") {
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
