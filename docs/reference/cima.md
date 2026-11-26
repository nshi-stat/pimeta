# Calculating Confidence Intervals

This function calculates confidence intervals.

## Usage

``` r
cima(
  y,
  se,
  v = NULL,
  alpha = 0.05,
  method = c("boot", "DL", "HK", "SJ", "KR", "APX", "PL", "BC"),
  B = 25000,
  parallel = FALSE,
  seed = NULL,
  maxit1 = 1e+05,
  eps = 10^(-10),
  lower = 0,
  upper = 1000,
  maxit2 = 1000,
  tol = .Machine$double.eps^0.25,
  rnd = NULL,
  maxiter = 100,
  sampling = c("table", "exact"),
  nodes = 400
)
```

## Arguments

- y:

  the effect size estimates vector

- se:

  the within studies standard errors vector

- v:

  the within studies variance estimates vector

- alpha:

  the alpha level of the prediction interval

- method:

  the calculation method for the confidence interval (default = "boot").

  - `boot`: A parametric bootstrap confidence interval (Nagashima et
    al., 2019).

  - `DL`: A Wald-type t-distribution confidence interval (the
    DerSimonian & Laird estimator for \\(\tau^2\\) with an approximate
    variance estimator for the average effect,
    \\((\sum{\hat{w}\_i})^{-1}\\), \\(df=K-1\\)).

  - `HK`: A Wald-type t-distribution confidence interval (the REML
    estimator for \\(\tau^2\\) with the Hartung (1999)'s variance
    estimator \[the Hartung and Knapp (2001)'s estimator\] for the
    average effect, \\(df=K-1\\)).

  - `SJ`: A Wald-type t-distribution confidence interval (the REML
    estimator for \\(\tau^2\\) with the Sidik and Jonkman (2006)'s bias
    corrected SE estimator for the average effect, \\(df=K-1\\)).

  - `KR`: Partlett–Riley (2017) confidence interval / (the REML
    estimator for \\(\tau^2\\) with the Kenward and Roger (1997)'s
    approach for the average effect, \\(df=\nu\\)).

  - `APX`: A Wald-type t-distribution confidence interval / (the REML
    estimator for \\(\tau^2\\) with an approximate variance estimator
    for the average effect, \\(df=K-1\\)).

  - `PL`: Profile likelihood confidence interval (Hardy & Thompson,
    1996).

  - `BC`: Profile likelihood confidence interval with Bartlett-type
    correction (Noma, 2011).

- B:

  the number of bootstrap samples

- parallel:

  the number of threads used in parallel computing, or FALSE that means
  single threading

- seed:

  set the value of random seed

- maxit1:

  the maximum number of iteration for the exact distribution function of
  \\(Q\\)

- eps:

  the desired level of accuracy for the exact distribution function of
  \\(Q\\)

- lower:

  the lower limit of random numbers of \\(\tau^2\\)

- upper:

  the upper limit of random numbers of \\(\tau^2\\) (see
  [`pima`](https://nshi-stat.github.io/pimeta/reference/pima.md))

- maxit2:

  the maximum number of iteration for numerical inversions

- tol:

  the desired level of accuracy for numerical inversions (see
  [`pima`](https://nshi-stat.github.io/pimeta/reference/pima.md))

- rnd:

  a vector of random numbers from the exact distribution of \\(\tau^2\\)

- maxiter:

  the maximum number of iteration for REML estimation

- sampling:

  the sampling method of \\(\tau^2\\) in the parametric bootstrap
  (`"table"` \[default\] or `"exact"`; see
  [`pima`](https://nshi-stat.github.io/pimeta/reference/pima.md))

- nodes:

  the number of nodes of the quantile table (see
  [`pima`](https://nshi-stat.github.io/pimeta/reference/pima.md))

## Value

- `K`: the number of studies.

- `muhat`: the average treatment effect estimate \\(\hat{\mu}\\).

- `lci`, `uci`: the lower and upper confidence limits \\(\hat{\mu}\_l\\)
  and \\(\hat{\mu}\_u\\).

- `tau2h`: the estimate for \\(\tau^2\\).

- `i2h`: the estimate for \\(I^2\\).

- `nuc`: degrees of freedom for the confidence interval.

- `vmuhat`: the variance estimate for \\(\hat{\mu}\\).

## Details

Excellent reviews of heterogeneity variance estimation have been
published (e.g., Veroniki, et al., 2019).

## References

Veroniki, A. A., Jackson, D., Bender, R., Kuss, O., Langan, D., Higgins,
J. P. T., Knapp, G., and Salanti, G. (2019). Methods to calculate
uncertainty in the estimated overall effect size from a random-effects
meta-analysis. *Res Synth Methods.* **10**(1): 23-43.
[doi:10.1002/jrsm.1319](https://doi.org/10.1002/jrsm.1319) .

Nagashima, K., Noma, H., and Furukawa, T. A. (2019). Prediction
intervals for random-effects meta-analysis: a confidence distribution
approach. *Stat Methods Med Res*. **28**(6): 1689-1702.
[doi:10.1177/0962280218773520](https://doi.org/10.1177/0962280218773520)
.

Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009). A
re-evaluation of random-effects meta-analysis. *J R Stat Soc Ser A Stat
Soc.* **172**(1): 137-159.
[doi:10.1111/j.1467-985X.2008.00552.x](https://doi.org/10.1111/j.1467-985X.2008.00552.x)

Partlett, C, and Riley, R. D. (2017). Random effects meta-analysis:
Coverage performance of 95 confidence and prediction intervals following
REML estimation. *Stat Med.* **36**(2): 301-317.
[doi:10.1002/sim.7140](https://doi.org/10.1002/sim.7140)

Hartung, J., and Knapp, G. (2001). On tests of the overall treatment
effect in meta-analysis with normally distributed responses. *Stat Med.*
**20**(12): 1771-1782.
[doi:10.1002/sim.791](https://doi.org/10.1002/sim.791)

Sidik, K., and Jonkman, J. N. (2006). Robust variance estimation for
random effects meta-analysis. *Comput Stat Data Anal.* **50**(12):
3681-3701.
[doi:10.1016/j.csda.2005.07.019](https://doi.org/10.1016/j.csda.2005.07.019)

Noma H. (2011) Confidence intervals for a random-effects meta-analysis
based on Bartlett-type corrections. *Stat Med.* **30**(28): 3304-3312.
[doi:10.1002/sim.4350](https://doi.org/10.1002/sim.4350)

## See also

[`pima`](https://nshi-stat.github.io/pimeta/reference/pima.md)

## Examples

``` r
data(sbp, package = "pimeta")
set.seed(20161102)

# Nagashima-Noma-Furukawa confidence interval
pimeta::cima(sbp$y, sbp$sigmak, seed = 3141592)
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A parametric bootstrap confidence interval
#>  Heterogeneity variance: DerSimonian-Laird
#>  Variance for average treatment effect: Hartung
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3341 [-0.5673, -0.0985]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0282
#>  I-squared: 70.5%
#> 

# A Wald-type t-distribution confidence interval
# An approximate variance estimator & DerSimonian-Laird estimator for tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "DL")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3341 [-0.5068, -0.1613]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0282
#>  I-squared: 70.5%
#> 

# A Wald-type t-distribution confidence interval
# The Hartung variance estimator & REML estimator for tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "HK")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A Wald-type t-distribution confidence interval
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: Hartung-Knapp
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5761, -0.0814]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared: 85.5%
#> 

# A Wald-type t-distribution confidence interval
# The Sidik-Jonkman variance estimator & REML estimator for tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "SJ")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A Wald-type t-distribution confidence interval
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: bias corrected Sidik-Jonkman
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5625, -0.0950]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared: 85.5%
#> 

# A Wald-type t-distribution confidence interval
# The Kenward-Roger approach & REML estimator for tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "KR")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A Wald-type t-distribution confidence interval
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: Kenward-Roger
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5815, -0.0759]
#>  d.f.: 6.9508
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared: 85.5%
#> 

# A Wald-type t-distribution confidence interval
# An approximate variance estimator & REML estimator for tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "APX")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5646, -0.0929]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared: 85.5%
#> 

# Profile likelihood confidence interval
# Maximum likelihood estimators of variance for the average effect & tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "PL")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A profile likelihood confidence interval
#>  Heterogeneity variance: ML
#>  Variance for average treatment effect: ML
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3301 [-0.5435, -0.1089]
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0547
#>  I-squared: 82.2%
#> 

# Profile likelihood confidence interval with a Bartlett-type correction
# Maximum likelihood estimators of variance for the average effect & tau^2
pimeta::cima(sbp$y, sbp$sigmak, method = "BC")
#> 
#> Confidence Interval for Random-Effects Meta-Analysis
#> 
#> A profile likelihood confidence interval with a Bartlett-type correction
#>  Heterogeneity variance: ML
#>  Variance for average treatment effect: ML
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3301 [-0.5733, -0.0781]
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0547
#>  I-squared: 82.2%
#> 
```
