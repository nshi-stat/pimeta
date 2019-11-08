# Prediction interval based on Bartlett type corrections
#
# Returns a confidence interval for \eqn{\hat{\mu}} based on
# Bartlett type corrections (Noma, 2011).
# 
# @name pima_wl
# @rdname pima_wl
# @param y the effect size estimates vector
# @param sigma the within studies standard errors vector
# @param alpha the alpha level of the prediction interval
# @return
# \itemize{
# \item \code{muhat}: the average treatment effect estimate \eqn{\hat{\mu}}.
# \item \code{lpi}, \code{lpi}: the lower and upper confidence limits
#                               \eqn{\hat{\mu}_l} and \eqn{\hat{\mu}_u}.
# \item \code{tau2h}: the estimate for \eqn{\tau^2}.
# }
# @references
# Wang, C-C and Lee, W-C. (2019).
# A simple method to estimate prediction intervals and predictive
# distributions.
# \emph{Res Syn Meth.}
# \strong{30}(28): 3304-3312.
# \url{https://doi.org/10.1002/sim.4350}.
# @seealso
# \code{\link[=pima]{pima}}.
# @examples
# data(sbp, package = "pimeta")
# pimeta::pima_wl(sbp$y, sbp$sigmak)
# @export
pima_wl <- function(y, sigma, alpha = 0.05) {
  
  # initial check
  util_check_num(y)
  util_check_nonneg(sigma)
  util_check_inrange(alpha, 0.0, 1.0)
  
  if (length(sigma) != length(y)) {
    stop("'y' and 'sigma' should have the same length.")
  }
  
  # estimation
  k <- length(y)
  tau2h <- tau2h_dl(y = y, se = sigma)$tau2h
  w <- (sigma^2 + tau2h)^-1
  muhat <- sum(y*w) / sum(w)
  vmuhat <- 1/sum(w)
  lci <- muhat - stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  uci <- muhat + stats::qt(1 - alpha*0.5, k - 1)*sqrt(vmuhat)
  
  new_mean <- y
  sd <- sigma

  ########################################
  # codes from Wang and Lee (2019).
  
  num <- length(new_mean)
  theta0 <- mean(new_mean)
  variance <- sd^2
  
  theta <- sum(new_mean*sd^(-2))/sum(sd^(-2))
  thetav <- rep(theta, num)
  Q <- sum((new_mean-thetav)^2*sd^(-2))
  v0 <- (Q-(num-1))/(sum(sd^(-2))-sum(sd^(-4))/sum(sd^(-2)))
  if (v0 < 0) {
    v0 <- 0
  }
  tau2 <- v0
  thetar <- sum(new_mean/(sd^2+tau2))/sum((sd^2+tau2)^(-1))
  post_mean <- rep(0, num)
  post_var <- rep(0, num)
  for (i in 1:num) {
    if (v0 > 0) {
      post_var[i] <- 1/(1/v0 + 1/variance[i])
      post_mean[i] <- thetar + (new_mean[i]-thetar)*(tau2/(tau2+variance[i]))^0.5
    }
    else {
      post_var[i] <- 0
      post_mean[i] <- thetar
    }
  }
  mean <- thetar
  post_sd_adjusted <- sqrt(post_var)
  post_mean_sorted <- sort(post_mean)
  cdfi <- function(x) {
    a <- floor(x*(num + 1))
    b <- a+1
    if (a == 0) {
      mu <- mean
      sigma <- (post_mean_sorted[1] - mu)/qnorm(1/(num + 1))
    } else if (a == num) {
      mu <- mean
      sigma <- (post_mean_sorted[num] - mu)/qnorm(num/(num + 1))
    }
    else {
      sigma <- (post_mean_sorted[b] - post_mean_sorted[a])/(qnorm(b/(num + 1)) - qnorm(a/(num + 1)))
      mu <- post_mean_sorted[a]- qnorm(a/(num + 1))*sigma
    }
    return(qnorm(x, mean = mu, sd = sigma))
  }
  if (v0 > 0) {
    PI.lower <- cdfi(alpha/2)
    PI.upper <- cdfi(1 - alpha/2)
  }
  else {
    PI.lower <- mean
    PI.upper <- mean
  }
  
  ########################################

  res <- list(muhat = thetar, lpi = PI.lower, upi = PI.upper, lci = lci, uci = uci,
              tau2h = tau2, vmuhat = vmuhat, nup = NA, nuc = k - 1,
              method = "WL", y = y, se = sigma, alpha = alpha)
  
  return(res)
  
}
