#' Calculating Confidence Intervals
#' 
#' This function calculates confidence intervals.
#'
#' @name cima
#' @rdname cima
#' @param y the effect size estimates vector
#' @param se the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @param method the calculation method for the pretiction interval (default = "boot").
#' \itemize{
#' \item \code{boot}: A parametric bootstrap confidenceS interval
#'            (Nagashima et al., 2018).
#' \item \code{DL}: A Wald-type confidence interval
#'            (the DerSimonian & Laird estimator for \eqn{\tau^2} with
#'            a standard SE estimator for the average effect,
#'            \eqn{(1/\sum{\hat{w}_i})^{-1}}).
#' \item \code{HK}: A Wald-type t-distribution confidence interval
#'            (the REML estimator for \eqn{\tau^2} with
#'            the Hartung and Knapp (2001)'s SE estimator for the average effect).
#' \item \code{SJ}: A Wald-type t-distribution confidence interval
#'            (the REML estimator for \eqn{\tau^2} with
#'            the Sidik and Jonkman (2006)'s bias coreccted SE estimator
#'            for the average effect).
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
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Nagashima, K., Noma, H., and Furukawa, T. A. (2018).
#' Prediction intervals for random-effects meta-analysis:
#' a confidence distribution approach.
#' \emph{Stat Methods Med Res}.
#' \emph{In press}.
#' \url{https://doi.org/10.1177/0962280218773520}.
#' @seealso
#' \code{\link[=pima_boot]{pima_boot()}}
#' @examples
#' data(sbp, package = "pimeta")
#' set.seed(20161102)
#' \donttest{pimeta::pima(sbp$y, sbp$sigmak, B = 25000)}
#' @export
cima <- function(y, se, alpha = 0.05, method = c("boot", "DL", "HK", "SJ"),
                 B = 25000, maxit1 = 100000, eps = 10^(-10), lower = 0, upper = 1000,
                 maxit2 = 1000, tol = .Machine$double.eps^0.25, rnd = NULL,
                 maxiter = 100) {
  
  # initial check
  lstm <- c("boot", "DL", "HK", "SJ")
  method <- match.arg(method)

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
  } else if (alpha < 0.0 | alpha > 1.0) {
    stop("'alpha' should be 0.0 < alpha < 1.0.")
  } else if (B < 1) {
    stop("'B' should be grater than 1.")
  }
  
  if (method == "boot") {
    if (B < 1000) {
      warning("'B' > 1000 is recommended.")
    }
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
  } else if (method == "DL") {
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
  }
  res <- append(res, list(i2h = i2h(se, res$tau2h)))
  class(res) <- "cima"
  
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
#' @method print cima
print.cima <- function(x, digits = 4, ...) {
  
  cat("\nPrediction Interval for Random-Effects Meta-Analysis\n\n")
  
  if (x$method == "boot") {
    cat(paste0("A parametric bootstrap confidence interval\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: Hartung\n\n"))
  } else if (x$method == "DL") {
    cat(paste0("A Wald-type interval\n",
               " Heterogeneity variance: DerSimonian-Laird\n",
               " SE for average treatment effect: standard\n\n"))
  } else if (x$method == "HK") {
    cat(paste0("A Wald-type t-distribution confidence interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Hartung-Knapp\n\n"))
  } else if (x$method == "SJ") {
    cat(paste0("A Wald-type t-distribution confidence interval\n",
               " Heterogeneity variance: REML\n",
               " SE for average treatment effect: Sidik-Jonkman\n\n"))
  }
  
  cat(paste0("No. of studies: ", length(x$y), "\n\n"))

  cat(paste0("Average treatment effect [", (1 - x$alpha)*100, "%CI]:\n"))
  cat(paste0(" ", format(round(x$muhat, digits), nsmall = digits), " [",
             format(round(x$lci, digits), nsmall = digits), ", ",
             format(round(x$uci, digits), nsmall = digits), "]\n\n"))
  
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
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method plot cima
plot.cima <- function(x, y = NULL, title = "Forest plot", base_size = 16,
                      base_family = "", ...) {
  
  idodr <- lcl <- limits <- lx <- shape <- size <- ucl <- ymax <- ymin <- NULL
  
  k <- length(x$y)
  id <- c(
    paste0("  ", 1:k),
    paste0("  95%CI (I^2 = ", sprintf("%4.1f", x$i2h), "%)")
  )
  df1 <- data.frame(
    id = id,
    idodr = c((k + 2):3, 1),
    y = c(x$y, NA),
    lcl = c(x$y + stats::qnorm(0.025)*x$se, NA),
    ucl = c(x$y + stats::qnorm(0.975)*x$se, NA),
    shape = c(rep(15, k), 18),
    swidth = c(rep(1, k), 3)
  )
  df1 <- data.frame(
    df1,
    size = c(1/x$se, 1),
    limits = c(paste0(sprintf("%3.2 f", df1$y[1:k]), " (",
                      sprintf("%3.2 f", df1$lcl[1:k]), ", ",
                      sprintf("%3.2 f", df1$ucl[1:k]), ")"),
               paste0(sprintf("%3.2 f", x$muhat), " (",
                      sprintf("%3.2 f", x$lci), ", ",
                      sprintf("%3.2 f", x$uci), ")")
    ),
    lx = rep(max(df1$ucl, na.rm = TRUE) + 1.7, k + 1)
  )
  df2 <- data.frame(id = "Study", idodr = k + 3, y = NA, lcl = NA, ucl = NA,
                    shape = NA, swidth = NA, size = NA, limits = NA, lx = NA)
  df3 <- data.frame(id = "Overall", idodr = 2, y = NA, lcl = NA, ucl = NA,
                    shape = NA, swidth = NA, size = NA, limits = NA, lx = NA)
  df4 <- data.frame(x = c(x$lci, x$muhat, x$uci), ymax = c(1, 1 + 0.25, 1),
                    ymin = c(1, 1 - 0.25, 1), y = c(1, 1, 1))
  df1 <- rbind(df3, df2, df1)
  
  ggplot <- ggplot2::ggplot
  aes <- ggplot2::aes
  geom_errorbarh <- ggplot2::geom_errorbarh
  geom_point <- ggplot2::geom_point
  geom_ribbon <- ggplot2::geom_ribbon
  geom_vline <- ggplot2::geom_vline
  geom_text <- ggplot2::geom_text
  scale_y_discrete <- ggplot2::scale_y_discrete
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
  
  suppressWarnings(print(
    p <- ggplot(df1, aes(x = y, y = stats::reorder(id, idodr))) +
      geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0, size = 1) +
      geom_point(aes(size = size, shape = shape), fill = "black", show.legend = FALSE) +
      geom_ribbon(data = df4, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 1,
                  colour = "black", fill = "black") +
      geom_vline(xintercept = x$muhat, lty = 2) +
      geom_vline(xintercept = 0, lty = 1) +
      geom_text(aes(label = limits, x = lx), size = base_size*0.282, hjust = 1) +
      scale_y_discrete() +
      scale_x_continuous() +
      scale_shape_identity() +
      ylab(NULL) +
      xlab("  ") +
      ggtitle(title) +
      theme_classic(base_size = base_size, base_family = "") +
      theme(axis.text.y = element_text(hjust = 0), axis.ticks.y = element_blank()) +
      theme(axis.line.x = element_line(), axis.line.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = rel(0.8)))
  ))
  
}