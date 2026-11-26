# Helpers for the unit tests.
#
# The reference implementations below are written independently of the
# package code, directly from the formulas in the original papers (see also
# Veroniki et al., 2016, Res Synth Methods 7: 55-79). The iterative
# estimators are obtained by maximising the (restricted) likelihood or by
# solving their estimating equations on a grid, not by the fixed-point
# iterations used in the package.

get_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "pimeta", envir = e)
  d <- e[[name]]
  list(y = d$y, se = if (is.null(d$sigmak)) d$se else d$sigmak)
}

datasets <- function() {
  lapply(c(sbp = "sbp", hyp = "hyp", pain = "pain", setshift = "setshift"),
         get_data)
}

# a homogeneous data set (all the estimators of tau^2 are 0)
homog <- list(y = c(0.10, 0.12, 0.09, 0.11, 0.10),
              se = c(0.20, 0.25, 0.30, 0.20, 0.22))

# |object - expected| <= tol * (1 + |expected|)
expect_close <- function(object, expected, tol = 1e-4) {
  testthat::expect_lte(abs(object - expected), tol*(1 + abs(expected)))
}

# generalised Q statistic
ref_q <- function(y, v, tau2 = 0) {
  w <- 1/(v + tau2)
  mu <- sum(w*y)/sum(w)
  sum(w*(y - mu)^2)
}

# a grid on [0, upper scale] used to locate roots and maxima
ref_grid <- function(y, v) {
  c(0, (stats::var(y) + max(v))*10^seq(-8, 4, length.out = 2000))
}

# the root of f at the last sign change from + to - on the grid
# (estimating equations; 0 if f <= 0 on the grid)
ref_root <- function(f, grid) {
  fv <- vapply(grid, f, 0.0)
  i <- which(fv[-length(fv)] > 0 & fv[-1] <= 0)
  if (length(i) == 0) return(0)
  i <- max(i)
  stats::uniroot(f, grid[c(i, i + 1)], tol = 1e-12)$root
}

# the maximiser of f on the grid, refined by optimize()
ref_max <- function(f, grid) {
  fv <- vapply(grid, f, 0.0)
  i <- which.max(fv)
  if (i == 1) return(0)
  stats::optimize(f, grid[c(i - 1, min(i + 1, length(grid)))],
                  maximum = TRUE, tol = 1e-12)$maximum
}

ref_tau2 <- function(y, se, method) {
  v <- se^2
  w0 <- 1/v
  k <- length(y)
  grid <- ref_grid(y, v)
  resid <- function(t) {
    w <- 1/(v + t)
    list(w = w, r2 = (y - sum(w*y)/sum(w))^2)
  }
  switch(method,
    DL = max(0, (ref_q(y, v) - (k - 1))/(sum(w0) - sum(w0^2)/sum(w0))),
    VC = max(0, stats::var(y) - mean(v)),
    HS = max(0, (ref_q(y, v) - k)/sum(w0)),
    HM = {
      q <- ref_q(y, v)
      q^2/((sum(w0) - sum(w0^2)/sum(w0))*(2*(k - 1) + q))
    },
    SJ = {
      t0 <- mean((y - mean(y))^2)
      t0*ref_q(y, v, t0)/(k - 1)
    },
    SJ2 = {
      t0 <- max(0.01, max(0, stats::var(y) - mean(v)))
      t0*ref_q(y, v, t0)/(k - 1)
    },
    PM = ref_root(function(t) ref_q(y, v, t) - (k - 1), grid),
    # ML and REML: maximisers of the (restricted) log-likelihood
    ML = ref_max(function(t) {
      a <- resid(t); -0.5*sum(log(v + t)) - 0.5*sum(a$w*a$r2)
    }, grid),
    REML = ref_max(function(t) {
      a <- resid(t); -0.5*sum(log(v + t)) - 0.5*sum(a$w*a$r2) - 0.5*log(sum(a$w))
    }, grid),
    AREML = ref_root(function(t) {
      a <- resid(t); k/(k - 1)*sum(a$w^2*a$r2) - sum(a$w)
    }, grid),
    EB = ref_root(function(t) {
      a <- resid(t); k/(k - 1)*sum(a$w*a$r2) - k
    }, grid)
  )
}

# log-likelihood and profile log-likelihood of mu
ref_ll <- function(mu, tau2, y, v) {
  -0.5*sum(log(v + tau2)) - 0.5*sum((y - mu)^2/(v + tau2))
}
ref_pll <- function(mu, y, v) {
  grid <- c(0, 10*(max((y - mu)^2) + max(v))*10^seq(-10, 0, length.out = 1000))
  f <- function(t) ref_ll(mu, t, y, v)
  f(ref_max(f, grid))
}
