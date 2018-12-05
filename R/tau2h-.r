#' Calculating Heterogeneity Variance
#' 
#' Returns a heterogeneity variance estimate and
#' its confidence interval.
#' 
#' @name tau2h
#' @rdname tau2h
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @param maxiter the maximum number of iterations
#' @param method the calculation method for heterogeneity
#' variance (default = "DL").
#' \itemize{
#' \item \code{DL}: DerSimonian & Laird (1986).
#' \item \code{VC}: Hedges (1983).
#' \item \code{PM}: Paule & Mandel (1982).
#' \item \code{HM}: Hedges (1983).
#' \item \code{HS}: Hedges (1983).
#' \item \code{ML}: Hedges (1983).
#' \item \code{REML}: Hedges (1983).
#' \item \code{AREML}: Hedges (1983).
#' \item \code{SJ}: Hedges (1983).
#' \item \code{SJ2}: Hedges (1983).
#' \item \code{EB}: Hedges (1983).
#' \item \code{BM}: Hedges (1983).
#' }
#' @return
#' \itemize{
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' \item \code{lci}, \code{uci}: the lower and upper confidence limits \eqn{\hat{\tau}^2_l} and \eqn{\hat{\tau}^2_u}.
#' }
#' @references
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
#' Mandel, R. C., and Mandel, K. H. (1982).
#' Consensus values and weighting factors. 
#' \emph{J Res Natl Inst Stand Techno.}
#' \strong{87}(5): 377-385. 
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::tau2h(sbp$y, sbp$sigmak, method = "DL")
#' @export
tau2h <- function(y, se, maxiter = 100, method = c("DL", "VC")) {
  
  ## .. need more more strictry check.
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (min(se) < 0.0) {
    stop("'se' should be positive.")
  }
  
  # tbd
  
}
