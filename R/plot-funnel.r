#' Funnel Plot
#' 
#' Function for funnel plot of `pima` or `cima` objects.
#'
#' @param x `pima` or `cima` object to plot
#' @param title graph title
#' @param base_size base font size
#' @param base_family base font family
#' @param digits a value for digits specifies the minimum number
#'               of significant digits to be printed in values.
#' @param trans transformation for logarithmic scale outcomes
#'              (\code{"identity"} [default] or \code{"exp"}).
#' @examples
#' \donttest{
#' data(sbp, package = "pimeta")
#' piex <- pimeta::pima(sbp$y, sbp$sigmak, method = "HTS")
#' cairo_pdf("forestplot.pdf", width = 5, height = 5, family = "Arial")
#' funnel(piex, digits = 2, base_size = 10)
#' dev.off()
#' }
#' @export
#' @method plot pima
funnel <- function(x, title = "Funnel plot", base_size = 16,
                   base_family = "", digits = 3, trans = c("identity", "exp")) {
  
  lstt <- c("identity", "exp")
  trans <- match.arg(trans)
  if (any(class(x) != "pima") | any(class(x) != "cima")) {
    stop("Argument 'x' must be an object of class \"pima\" or \"cima\".")
  } else if (!is.element(trans, lstt)) {
    stop("Unknown 'trans' specified.")
  }
  
  ftrans <- function(x) {
    if (trans == "identity") {
      return(x)
    } else if (trans == "exp") {
      return(exp(x))
    }
  }
  
  dat <- data.frame(y = x$y, se = x$se)

  
}
