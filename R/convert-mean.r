#' Converting means and standard deviations
#'
#' Converting the means and standard deviations in the experimental
#' and control groups to the effect size estimates (mean differences)
#' and the within studies standard errors vector
#'
#' @details The effect size is the mean difference
#' \eqn{y = m1 - m2}. The within study standard error is
#' \eqn{\sqrt{s1^2/n1 + s2^2/n2}} (\code{pooled = FALSE}, default) or
#' \eqn{s_p \sqrt{1/n1 + 1/n2}} with the pooled standard deviation
#' \eqn{s_p^2 = ((n1-1) s1^2 + (n2-1) s2^2)/(n1+n2-2)}
#' (\code{pooled = TRUE}).
#'
#' @name convert_mean
#' @rdname convert_mean
#' @param n1 a vector of the number of patients in the experimental group
#' @param m1 a vector of the means in the experimental group
#' @param s1 a vector of the standard deviations in the experimental group
#' @param n2 a vector of the number of patients in the control group
#' @param m2 a vector of the means in the control group
#' @param s2 a vector of the standard deviations in the control group
#' @param pooled logical; if \code{TRUE}, the pooled variance is used.
#'   The default is \code{FALSE}.
#' @return A \code{data.frame} of study data.
#' \itemize{
#' \item \code{y}: a numeric vector of the effect size estimates
#'   (mean differences).
#' \item \code{se}: a numeric vector of the within studies standard errors.
#' }
#' @seealso \code{\link{convert_bin}}, \code{\link{clbp}}, \code{\link{nfld}}
#' @examples
#' data(clbp, package = "pimeta")
#' dat <- pimeta::convert_mean(clbp$n1, clbp$m1, clbp$s1,
#'                             clbp$n2, clbp$m2, clbp$s2)
#' head(dat)
#' \donttest{pimeta::pima(dat$y, dat$se, method = "HTS")}
#' @export
convert_mean <- function(n1, m1, s1, n2, m2, s2, pooled = FALSE) {

  # initial check
  util_check_nonneg(n1)
  util_check_num(m1)
  util_check_nonneg(s1)
  util_check_nonneg(n2)
  util_check_num(m2)
  util_check_nonneg(s2)

  if (!is.logical(pooled) || length(pooled) != 1 || is.na(pooled)) {
    stop("'pooled' should be TRUE or FALSE.")
  }
  if (length(n1) != length(m1) || length(m1) != length(s1) || length(s1) != length(n2) ||
      length(n2) != length(m2) || length(m2) != length(s2)) {
    stop("'n1', 'm1', 's1', 'n2', 'm2', and 's2' should have the same length.")
  }
  if (any(n1 <= 0) || any(n2 <= 0)) {
    stop("'n1' and 'n2' should be positive.")
  }
  if (pooled && any(n1 + n2 <= 2)) {
    stop("'n1 + n2' should be greater than 2 when 'pooled = TRUE'.")
  }

  res <- NULL
  res$y <- m1 - m2
  if (pooled) {
    sp2 <- ((n1 - 1)*s1^2 + (n2 - 1)*s2^2)/(n1 + n2 - 2)
    res$se <- sqrt(sp2*(1/n1 + 1/n2))
  } else {
    res$se <- sqrt(s1^2/n1 + s2^2/n2)
  }

  res <- data.frame(res)

  return(res)

}
