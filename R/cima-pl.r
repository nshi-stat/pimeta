#' Brockwell--Gordon confidence interval
#'
#' Returns the Brockwell--Gordon confidence interval
#' for \eqn{\hat{\mu}} (Brockwell & Gordon, 2007).
#' 
#' @name cima_pl
#' @rdname cima_pl
#' @aliases htsdl
#' @param y the effect size estimates vector
#' @param sigma the within studies standard errors vector
#' @param alpha the alpha level of the prediction interval
#' @return
#' \itemize{
#' \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
#' \item \code{lci}, \code{lci}: the lower and upper confidence limits
#'                               \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
#' \item \code{tau2h}: the estimate for \eqn{\tau^2}.
#' }
#' @references
#' Hardy, R. J., and Thompson, S. G. (1996).
#' A likelihood approach to meta-analysis with random effects.
#' \emph{Stat Med.}
#' \strong{15}(6): 619-629.
#' \url{https://doi.org/10.1002/(SICI)1097-0258(19960330)15:6<619::AID-SIM188>3.0.CO;2-A}.
#' @seealso
#' \code{\link[=cima]{cima}}.
#' @examples
#' data(sbp, package = "pimeta")
#' pimeta::cima_pl(sbp$y, sbp$sigmak)
#' @export
cima_pl <- function(y, se, alpha = 0.05) {

  # tbd
  if (0) {
    
  ll <- function(mu, tau2, y, se, k) {
    -0.5*log(2.0*pi) - 0.5*sum(log(se^2 + tau2)) - 0.5*sum((y - mu)^2/(se^2 + tau2))
  }
  
  k <- length(y)
  tau2h <- tau2h_ml(y = y, se = se)$tau2h
  w <- (se^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  
  llp <- function(mu, tau2, y, se, k) {
    
  }
  
  res <- try(
    uniroot(qgen, interval = c(muhat, upper), y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  tau2h <- res$root
  
  w <- (se^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  vmuhat <- 1/sum(w)
  lci <- muhat - stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  res <- list(muhat = muhat, lci = lci, uci = uci, tau2h = tau2h,
              method = "DLt", y = y, se = se, alpha = alpha)

  return(res)
  }
  
}

