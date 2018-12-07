<<<<<<< HEAD
#' Profile likelihood confidence interval
#'
#' Returns a profile likelihood confidence interval
#' for \eqn{\hat{\mu}} (Hardy & Thompson, 1996).
=======
#' Brockwell--Gordon confidence interval
#'
#' Returns the Brockwell--Gordon confidence interval
#' for \eqn{\hat{\mu}} (Brockwell & Gordon, 2007).
>>>>>>> 2efdcfa597cf8a71a55b49e1d047d17a1dd22c0f
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
<<<<<<< HEAD
=======

  # tbd
  if (0) {
    
  ll <- function(mu, tau2, y, se, k) {
    -0.5*log(2.0*pi) - 0.5*sum(log(se^2 + tau2)) - 0.5*sum((y - mu)^2/(se^2 + tau2))
  }
>>>>>>> 2efdcfa597cf8a71a55b49e1d047d17a1dd22c0f
  
  k <- length(y)
  tau2h <- tau2h_ml(y = y, se = se)$tau2h
  w <- (se^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  
<<<<<<< HEAD
  upper <- abs(muhat)*100
  while(peqn(upper, muhat, alpha, y, se) > 0) {
    upper <- upper*2
  }
  res <- try(
    uniroot(peqn, interval = c(muhat, upper), muhat = muhat,
            alpha = alpha, y = y, se = se),
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
  uci <- res$root
  
  lower <- -abs(muhat)*100
  while(peqn(lower, muhat, alpha, y, se) > 0) {
    lower <- lower*2
  }
  res <- try(
    uniroot(peqn, interval = c(lower, muhat), muhat = muhat,
            alpha = alpha, y = y, se = se),
=======
  llp <- function(mu, tau2, y, se, k) {
    
  }
  
  res <- try(
    uniroot(qgen, interval = c(muhat, upper), y = y, se = se),
>>>>>>> 2efdcfa597cf8a71a55b49e1d047d17a1dd22c0f
    silent = FALSE
  )
  if (class(res) == "try-error") {
    stop("Could not find a solution. Try another estimator.")
  }
<<<<<<< HEAD
  lci <- res$root
  
  res <- list(muhat = muhat, lci = lci, uci = uci, tau2h = tau2h,
              method = "PL", y = y, se = se, alpha = alpha)
  
  return(res)
  
}

peqn <- function(mu, muhat, alpha, y, se) {
  lpl(mu, y, se) - lpl(muhat, y, se) + 0.5*stats::qchisq(1.0 - alpha, 1)
}

lpl <- function(mu, y, se) {
  tau2mu <- tau2h_pl(mu, y, se)
  -0.5*sum(log(se^2 + tau2mu)) - 0.5*sum((y - mu)^2/(se^2 + tau2mu))
}

tau2h_pl <- function(mu, y, se, maxiter = 100) {
  
  tau2h <- tau2h_dl(y, se)$tau2h
  r <- 0
  
  # auto ajdustment for step length
  autoadj <- 0
  stepadj <- 0.5
  while(1) {
    wi <- (se^2 + tau2h)^-1
    ti <- sum(wi^2 * ((y - mu)^2 - se^2)) / sum(wi^2)
    if (ti <= 0) {
      tau2h <- 0.0
      break
    } else {
      if (abs(ti - tau2h)/(1.0 + tau2h) < 1e-5) {
        # converged
        break
      } else if (r == maxiter && autoadj == 1) {
        # not converged
        stop("The heterogeneity variance (tau^2) could not be calculated.")
        break
      } else if (r == maxiter && autoadj == 0) {
        # not converged but retry with an adjusted step length
        r <- 0
        autoadj <- 1
      } else if (autoadj == 1) {
        # with the adjustment
        r <- r + 1
        tau2h <- tau2h + (ti - tau2h)*stepadj
      } else {
        # without the adjustment
        r <- r + 1
        tau2h <- ti
      }
    }
  }
  
  return(tau2h)
  
=======
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
  
>>>>>>> 2efdcfa597cf8a71a55b49e1d047d17a1dd22c0f
}

