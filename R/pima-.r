#' Calculating Prediction Intervals
#' 
#' This function calculates prediction intervals using \code{pima_boot},
#' \code{pima_hts}, or \code{pima_htsreml} functions.
#'
#' @name pima
#' @rdname pima
#' @param y the effect size estimates vector
#' @param se the within studies standard error estimates vector
#' @param v the within studies variance estimates vector
#' @param alpha the alpha level of the prediction interval
#' @param method the calculation method for the pretiction interval (default = "boot").
#' \itemize{
#' \item \code{boot}: A parametric bootstrap prediction interval
#'            (Nagashima et al., 2018).
#' \item \code{HTS}: the Higgins--Thompson--Spiegelhalter (2009) prediction interval / 
#'            (the DerSimonian & Laird estimator for \eqn{\tau^2} with
#'            an approximate SE estimator for the average effect,
#'            \eqn{(1/\sum{\hat{w}_i})^{-1}}).
#' \item \code{HK}: Partlett--Riley (2017) prediction interval
#'            (the REML estimator for \eqn{\tau^2} with
#'            the Hartung and Knapp (2001)'s SE estimator for the average effect).
#' \item \code{SJ}: Partlett--Riley (2017) prediction interval /
#'            (the REML estimator for \eqn{\tau^2} with
#'            the Sidik and Jonkman (2006)'s bias coreccted SE estimator
#'            for the average effect).
#' \item \code{KR}: Partlett--Riley (2017) prediction interval /
#'            (the REML estimator for \eqn{\tau^2} with
#'            the Kenward and Roger (1997)'s approach
#'            for the average effect).
#' \item \code{APX}: Partlett--Riley (2017) prediction interval /
#'            (the REML estimator for \eqn{\tau^2} with
#'            an approximate SE estimator for the average effect).
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
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{uci}: the lower and upper confidence limits \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{lpi}, \code{upi}: the lower and upper prediction limits \eqn{\hat{c}_l} and \eqn{\hat{c}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' \item \code{vmuhat}: the variance estimate for \eqn{\hat{\mu}}.
#' \item \code{nup}: degrees of freedom for the prediction interval.
#' \item \code{nuc}: degrees of freedom for the confidence interval.
#' }
#' @references
#' Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009).
#' A re-evaluation of random-effects meta-analysis.
#' \emph{J R Stat Soc Ser A Stat Soc.}
#' \strong{172}(1): 137-159.
#' \url{https://doi.org/10.1111/j.1467-985X.2008.00552.x}
#'  
#' Partlett, C, and Riley, R. D. (2017).
#' Random effects meta-analysis: Coverage performance of 95%
#' confidence and prediction intervals following REML estimation.
#' \emph{Stat Med.}
#' \strong{36}(2): 301-317.
#' \url{https://doi.org/10.1002/sim.7140}
#' 
#' Nagashima, K., Noma, H., and Furukawa, T. A. (2018).
#' Prediction intervals for random-effects meta-analysis:
#' a confidence distribution approach.
#' \emph{Stat Methods Med Res}.
#' \emph{In press}.
#' \url{https://doi.org/10.1177/0962280218773520}.
#' @seealso
#' \code{\link[=pima_boot]{pima_boot}},
#' \code{\link[=pima_hts]{pima_hts}},
#' \code{\link[=pima_htsreml]{pima_htsreml}}.
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, B = 50000)}
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # A parametric bootstrap prediction interval
#' #  Heterogeneity variance: DerSimonian-Laird
#' #  SE for average treatment effect: Hartung
#' #
#' # No. of studies: 10
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3341 [-0.8715, 0.2144]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3341 [-0.5655, -0.0988]
#' # 
#' # Heterogeneity measure
#' #  tau2: 0.0282
#' #  I^2:  70.5%
#'  
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "HTS")}
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # Higgins-Thompson-Spiegelhalter prediction interval
#' #  Heterogeneity variance: DerSimonian-Laird
#' #  SE for average treatment effect: approximate
#' # 
#' # No. of studies: 10
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3341 [-0.7598, 0.0917]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3341 [-0.5068, -0.1613]
#' # 
#' # Heterogeneity measure
#' #  tau2: 0.0282
#' #  I^2:  70.5%
#' 
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "HK")}
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # Partlett-Riley prediction interval
#' #  Heterogeneity variance: REML
#' #  SE for average treatment effect: Hartung-Knapp
#' # 
#' # No. of studies: 10
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3287 [-0.9887, 0.3312]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3287 [-0.5761, -0.0814]
#' # 
#' # Heterogeneity measure
#' #  tau2: 0.0700
#' #  I^2:  85.5%
#' 
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, method = "SJ")}
#' # 
#' # Prediction Interval for Random-Effects Meta-Analysis
#' # 
#' # Partlett-Riley prediction interval
#' # Heterogeneity variance: REML
#' # SE for average treatment effect: Sidik-Jonkman
#' # 
#' # No. of studies: 10
#' # 
#' # Average treatment effect [95%PI]:
#' #  -0.3287 [-0.9835, 0.3261]
#' # 
#' # Average treatment effect [95%CI]:
#' #  -0.3287 [-0.5625, -0.0950]
#' # 
#' # Heterogeneity measure
#' #  tau2: 0.0700
#' #  I^2:  85.5%
#'
#' @export
pima <- function(y, se, v = NULL, alpha = 0.05,
                 method = c("boot", "HTS", "HK", "SJ", "KR", "CL", "APX"),
                 B = 25000, maxit1 = 100000, eps = 10^(-10), lower = 0, upper = 1000,
                 maxit2 = 1000, tol = .Machine$double.eps^0.25, rnd = NULL,
                 maxiter = 100) {
  
  # initial check
  lstm <- c("boot", "HTS", "HK", "SJ", "KR", "CL", "APX")
  method <- match.arg(method)
  
  if (missing(se) & missing(v)) {
    stop("Either 'se' or 'v' must be specified.")
  } else if (missing(se)) {
    se <- sqrt(v)
  }
  
  util_check_num(y)
  util_check_num(se)
  util_check_num(alpha)
  util_check_num(B)
  util_check_num(maxit1)
  util_check_num(eps)
  util_check_num(lower)
  util_check_num(upper)
  util_check_num(maxit2)
  util_check_num(tol)
  util_check_num(maxiter)
  util_check_nonneg(se)
  util_check_inrange(alpha, 0.0, 1.0)
  util_check_gt(B, 1)
  util_check_gt(maxit1, 1)
  util_check_gt(eps, 0)
  util_check_ge(lower, 0)
  util_check_gt(upper, 0)
  util_check_gt(maxit2, 1)
  util_check_gt(tol, 0)
  util_check_gt(maxiter, 1)
  
  if (length(se) != length(y)) {
    stop("'y' and 'se' should have the same length.")
  } else if (!is.element(method, lstm)) {
    stop("Unknown 'method' specified.")
  } else if (lower >= upper) {
    stop("'upper' should be greater than 'lower'.")
  }
  
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
  } else if (method == "KR") {
    res <- pima_htsreml(y       = y, 
                        sigma   = se, 
                        alpha   = alpha,
                        vartype = "KR",
                        maxiter = maxiter)
  } else if (method == "CL" | method == "APX") {
    res <- pima_htsreml(y       = y, 
                        sigma   = se, 
                        alpha   = alpha,
                        vartype = "APX",
                        maxiter = maxiter)
  }
  res <- append(res, list(i2h = i2h(se, res$tau2h)))
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
  
  nuc <- x$nuc
  nup <- x$nup
  
  cat("\nPrediction Interval for Random-Effects Meta-Analysis\n\n")
  
  if (x$method == "boot") {
    cat(paste0("A parametric bootstrap prediction and confidence intervals\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: Hartung\n\n"))
  } else if (x$method == "HTS") {
    cat(paste0("Higgins-Thompson-Spiegelhalter prediction interval\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: approximate\n\n"))
  } else if (x$method == "HK") {
    cat(paste0("Partlett-Riley prediction interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Hartung-Knapp\n\n"))
  } else if (x$method == "SJ") {
    cat(paste0("Partlett-Riley prediction interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: bias corrected Sidik-Jonkman\n\n"))
  } else if (x$method == "KR") {
    cat(paste0("Partlett-Riley prediction interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Kenward-Roger\n\n"))
    nup <- format(round(nup, digits))
    nuc <- format(round(nuc, digits))
  } else if (x$method == "CL" | x$method == "APX") {
    cat(paste0("A prediction interval with REML and standard SE\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: approximate\n\n"))
  }
  
  cat(paste0("No. of studies: ", length(x$y), "\n\n"))
  
  cat(paste0("Average treatment effect [", (1 - x$alpha)*100, "%PI]:\n"))
  cat(paste0(" ", format(round(x$muhat, digits), nsmall = digits), " [",
             format(round(x$lpi, digits), nsmall = digits), ", ",
             format(round(x$upi, digits), nsmall = digits), "]\n"))
  cat(paste0(" d.f.: ", nup, "\n\n"))
  
  cat(paste0("Average treatment effect [", (1 - x$alpha)*100, "%CI]:\n"))
  cat(paste0(" ", format(round(x$muhat, digits), nsmall = digits), " [",
             format(round(x$lci, digits), nsmall = digits), ", ",
             format(round(x$uci, digits), nsmall = digits), "]\n"))
  cat(paste0(" d.f.: ", nuc, "\n\n"))
  
  cat(paste0("Heterogeneity measure\n"))
  cat(paste0(" tau2: ", format(round(x$tau2, digits), nsmall = digits), "\n"))
  cat(paste0(" I^2:  ", format(round(x$i2h, 1), nsmall = 1), "%\n\n"))
  
  invisible(x)
  
}


#' Plot Results
#' 
#' A function for plotting of `pima` objects.
#'
#' @param x `pima` object to plot
#' @param y is not used
#' @param title graph title
#' @param base_size base font size
#' @param base_family base font family
#' @param digits a value for digits specifies the minimum number
#'               of significant digits to be printed in values.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method plot pima
plot.pima <- function(x, y = NULL, title = "Forest plot", base_size = 16,
                      base_family = "", digits = 3, studylabel = NULL, ...) {
  
  idodr <- lcl <- limits <- lx <- shape <- size <- ucl <- ymax <- ymin <- NULL
  
  k <- length(x$y)
  
  if (is.null(studylabel)) {
    studylabel <- 1:k
  } else {
    if (k != length(studylabel)) {
      stop("`studylabel` and the number of studies must have the same length.")
    }
  }
  
  id <- c(
    paste0("  ", studylabel),
    paste0("  95%CI (I^2 = ", format(round(x$i2h, 1), nsmall = 1), "%)"),
    paste0("  95%PI")
  )
  df1 <- data.frame(
    id = id,
    idodr = c((k + 3):4, 2:1),
    y = c(x$y, NA, NA),
    lcl = c(x$y + stats::qnorm(0.025)*x$se, NA, NA),
    ucl = c(x$y + stats::qnorm(0.975)*x$se, NA, NA),
    shape = c(rep(15, k), 18, 18),
    swidth = c(rep(1, k), 3, 3)
  )
  df1 <- data.frame(
    df1,
    size = c(1/x$se, 1, 1),
    limits = c(paste0(format(round(df1$y[1:k], digits), nsmall = digits), " (",
                      format(round(df1$lcl[1:k], digits), nsmall = digits), ", ",
                      format(round(df1$ucl[1:k], digits), nsmall = digits), ")"),
               paste0(format(round(x$muhat, digits), nsmall = digits), " (",
                      format(round(x$lci, digits), nsmall = digits), ", ",
                      format(round(x$uci, digits), nsmall = digits), ")"),
               paste0(format(round(x$muhat, digits), nsmall = digits), " (",
                      format(round(x$lpi, digits), nsmall = digits), ", ",
                      format(round(x$upi, digits), nsmall = digits), ")")
    )
  )
  df2 <- data.frame(id = "Study", idodr = k + 4, y = NA, lcl = NA, ucl = NA,
                    shape = NA, swidth = NA, size = NA, limits = NA)
  df3 <- data.frame(id = "Overall", idodr = 3, y = NA, lcl = NA, ucl = NA,
                    shape = NA, swidth = NA, size = NA, limits = NA)
  df4 <- data.frame(x = c(x$lci, x$muhat, x$uci), ymax = c(2, 2 + 0.25, 2),
                    ymin = c(2, 2 - 0.25, 2), y = c(2, 2, 2))
  df5 <- data.frame(x = c(x$lpi, x$muhat, x$upi), ymax = c(1, 1 + 0.25, 1),
                    ymin = c(1, 1 - 0.25, 1), y = c(1, 1, 1))
  df1 <- rbind(df3, df2, df1)
  
  ggplot <- ggplot2::ggplot
  aes <- ggplot2::aes
  geom_errorbarh <- ggplot2::geom_errorbarh
  geom_point <- ggplot2::geom_point
  geom_ribbon <- ggplot2::geom_ribbon
  geom_vline <- ggplot2::geom_vline
  scale_y_continuous <- ggplot2::scale_y_continuous
  scale_x_continuous <- ggplot2::scale_x_continuous
  scale_shape_identity <- ggplot2::scale_shape_identity
  ylab <- ggplot2::ylab
  xlab <- ggplot2::xlab
  ggtitle <- ggplot2::ggtitle
  theme_classic <- ggplot2::theme_classic
  theme <- ggplot2::theme
  element_text <- ggplot2::element_text
  element_line <- ggplot2::element_line
  element_blank <- ggplot2::element_blank
  element_blank <- ggplot2::element_blank
  rel <- ggplot2::rel
  labs <- ggplot2::labs
  sec_axis <- ggplot2::sec_axis
  
  y1labels <- df1[order(df1$idodr),]$id
  y2labels <- df1[order(df1$idodr),]$limits
  y2labels[is.na(y2labels)] <- ""
  
  suppressWarnings(
    print(
      p <- ggplot(df1, aes(x = y, y = idodr)) +
        geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0, size = 1) +
        geom_point(aes(size = size, shape = shape), fill = "black", show.legend = FALSE) +
        geom_ribbon(data = df4, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 1,
                    colour = "black", fill = "black") +
        geom_ribbon(data = df5, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 1,
                    colour = "black", fill = "black") +
        geom_vline(xintercept = x$muhat, lty = 2) +
        geom_vline(xintercept = 0, lty = 1) +
        scale_y_continuous(breaks = 1:length(y1labels), labels = y1labels,
                           sec.axis = sec_axis( ~ ., breaks = 1:length(y2labels), labels = y2labels)) +
        scale_x_continuous() +
        scale_shape_identity() +
        ylab(NULL) +
        xlab("  ") +
        labs(caption = parse(text = sprintf('hat(tau)^{2}=="%s"', format(round(x$tau2h, digits), 
                                                                         nsmall = digits)))) +
        ggtitle(title) +
        theme_classic(base_size = base_size, base_family = base_family) +
        theme(axis.text.y = element_text(hjust = 0), axis.ticks.y = element_blank()) +
        theme(axis.line.x = element_line(), axis.line.y = element_blank(),
              plot.title = element_text(hjust = 0.5, size = rel(0.8)))
    )
  )
  
}
