# Calculating Prediction Intervals

This function can estimate prediction intervals (PIs) as follows: A
parametric bootstrap PI based on confidence distribution (Nagashima et
al., 2019). A parametric bootstrap confidence interval is also
calculated based on the same sampling method for bootstrap PI. The
Higgins–Thompson–Spiegelhalter (2009) prediction interval. The
Partlett–Riley (2017) prediction intervals.

## Usage

``` r
pima(
  y,
  se,
  v = NULL,
  alpha = 0.05,
  method = c("boot", "HTS", "HK", "SJ", "KR", "CL", "APX"),
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

  the within studies standard error estimates vector

- v:

  the within studies variance estimates vector

- alpha:

  the alpha level of the prediction interval

- method:

  the calculation method for the prediction interval (default = "boot").

  - `boot`: A parametric bootstrap prediction interval (Nagashima et
    al., 2019).

  - `HTS`: the Higgins–Thompson–Spiegelhalter (2009) prediction interval
    / (the DerSimonian & Laird estimator for \\(\tau^2\\) with an
    approximate variance estimator for the average effect,
    \\((\sum{\hat{w}\_i})^{-1}\\), \\(df=K-2\\)).

  - `HK`: Partlett–Riley (2017) prediction interval (the REML estimator
    for \\(\tau^2\\) with the Hartung (1999)'s variance estimator \[the
    Hartung and Knapp (2001)'s estimator\] for the average effect,
    \\(df=K-2\\)).

  - `SJ`: Partlett–Riley (2017) prediction interval / (the REML
    estimator for \\(\tau^2\\) with the Sidik and Jonkman (2006)'s bias
    coreccted variance estimator for the average effect, \\(df=K-2\\)).

  - `KR`: Partlett–Riley (2017) prediction interval / (the REML
    estimator for \\(\tau^2\\) with the Kenward and Roger (1997)'s
    approach for the average effect, \\(df=\nu-1\\)).

  - `APX`: Partlett–Riley (2017) prediction interval / (the REML
    estimator for \\(\tau^2\\) with an approximate variance estimator
    for the average effect, \\(df=K-2\\)).

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

  the upper limit of random numbers of \\(\tau^2\\); the search is
  automatically extended to at least `1e8` times the mean of `sigma^2`
  (and `1e6` \\(Q\\) times). Random numbers beyond it are treated as
  \\(\tau^2 = \infty\\) (with a warning).

- maxit2:

  the maximum number of iteration for numerical inversions

- tol:

  the desired level of accuracy for numerical inversions, relative to
  the mean of `sigma^2` (the confidence distribution of \\(\tau^2\\) is
  computed on this standardised scale, so that the results are
  equivariant to the scale of `y` and `sigma`)

- rnd:

  a vector of random numbers from the exact distribution of \\(\tau^2\\)

- maxiter:

  the maximum number of iteration for REML estimation

- sampling:

  the sampling method of \\(\tau^2\\) from its confidence distribution
  in the parametric bootstrap (`method = "boot"`).

  - `table` (default): the inverse of the confidence distribution is
    computed exactly at `nodes - 1` probabilities and the bootstrap
    samples are obtained by monotone cubic (Fritsch–Carlson)
    interpolation of this quantile table. Probabilities outside the
    table are inverted exactly. The interpolation error is of the same
    order as the default accuracy of the numerical inversion (`tol`),
    while the number of numerical inversions is reduced from `B` to
    about `nodes`.

  - `exact`: the confidence distribution is inverted numerically for
    each bootstrap sample (the method of pimeta 1.1.3 and earlier).

  When `lower > 0`, `"exact"` is always used.

- nodes:

  the number of nodes of the quantile table (`sampling = "table"` only).

## Value

- `K`: the number of studies.

- `muhat`: the average treatment effect estimate \\(\hat{\mu}\\).

- `lci`, `uci`: the lower and upper confidence limits \\(\hat{\mu}\_l\\)
  and \\(\hat{\mu}\_u\\).

- `lpi`, `upi`: the lower and upper prediction limits \\(\hat{c}\_l\\)
  and \\(\hat{c}\_u\\).

- `tau2h`: the estimate for \\(\tau^2\\).

- `i2h`: the estimate for \\(I^2\\).

- `nup`: degrees of freedom for the prediction interval.

- `nuc`: degrees of freedom for the confidence interval.

- `vmuhat`: the variance estimate for \\(\hat{\mu}\\).

- `rnd`: the bootstrap samples of \\(\tau^2\\) (`method = "boot"` only).

## Details

The functions `bootPI`, `pima_boot`, `pima_hts`, `htsdl`,
`pima_htsreml`, `htsreml` are deprecated, and integrated to the `pima`
function.

## References

Higgins, J. P. T, Thompson, S. G., Spiegelhalter, D. J. (2009). A
re-evaluation of random-effects meta-analysis. *J R Stat Soc Ser A Stat
Soc.* **172**(1): 137-159.
<https://doi.org/10.1111/j.1467-985X.2008.00552.x>

Partlett, C, and Riley, R. D. (2017). Random effects meta-analysis:
Coverage performance of 95 confidence and prediction intervals following
REML estimation. *Stat Med.* **36**(2): 301-317.
<https://doi.org/10.1002/sim.7140>

Nagashima, K., Noma, H., and Furukawa, T. A. (2019). Prediction
intervals for random-effects meta-analysis: a confidence distribution
approach. *Stat Methods Med Res*. **28**(6): 1689-1702.
<https://doi.org/10.1177/0962280218773520>.

Hartung, J. (1999). An alternative method for meta-analysis. *Biom J.*
**41**(8): 901-916.
<https://doi.org/10.1002/(SICI)1521-4036(199912)41:8%3C901::AID-BIMJ901%3E3.0.CO;2-W>.

Hartung, J., and Knapp, G. (2001). On tests of the overall treatment
effect in meta-analysis with normally distributed responses. *Stat Med.*
**20**(12): 1771-1782. <https://doi.org/10.1002/sim.791>.

Sidik, K., and Jonkman, J. N. (2006). Robust variance estimation for
random effects meta-analysis. *Comput Stat Data Anal.* **50**(12):
3681-3701. <https://doi.org/10.1016/j.csda.2005.07.019>.

Kenward, M. G., and Roger, J. H. (1997). Small sample inference for
fixed effects from restricted maximum likelihood. *Biometrics.*
**53**(3): 983-997. <https://www.ncbi.nlm.nih.gov/pubmed/9333350>.

DerSimonian, R., and Laird, N. (1986). Meta-analysis in clinical trials.
*Control Clin Trials.* **7**(3): 177-188.

## See also

[`print.pima`](https://nshi-stat.github.io/pimeta/reference/print.pima.md),
[`plot.pima`](https://nshi-stat.github.io/pimeta/reference/plot.pima.md),
[`cima`](https://nshi-stat.github.io/pimeta/reference/cima.md).

## Examples

``` r
data(sbp, package = "pimeta")

# Nagashima-Noma-Furukawa prediction interval
# is sufficiently accurate when I^2 >= 10% and K >= 3
pimeta::pima(sbp$y, sbp$sigmak, seed = 3141592, parallel = 4)
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> A parametric bootstrap prediction and confidence intervals
#>  Heterogeneity variance: DerSimonian-Laird
#>  Variance for average treatment effect: Hartung (Hartung-Knapp)
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3341 [-0.8789, 0.2165]
#>  d.f.: 9
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3341 [-0.5673, -0.0985]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0282
#>  I-squared:  70.5%
#> 

# Higgins-Thompson-Spiegelhalter prediction interval and
# Partlett-Riley prediction intervals
# are accurate when I^2 > 30% and K > 25
pimeta::pima(sbp$y, sbp$sigmak, method = "HTS")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Higgins-Thompson-Spiegelhalter prediction and confidence intervals
#>  Heterogeneity variance: DerSimonian-Laird
#>  Variance for average treatment effect: approximate
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3341 [-0.7598, 0.0917]
#>  d.f.: 8
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3341 [-0.5068, -0.1613]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0282
#>  I-squared:  70.5%
#> 
pimeta::pima(sbp$y, sbp$sigmak, method = "HK")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Partlett-Riley prediction and confidence intervals
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: Hartung (Hartung-Knapp)
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3287 [-0.9887, 0.3312]
#>  d.f.: 8
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5761, -0.0814]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared:  85.5%
#> 
pimeta::pima(sbp$y, sbp$sigmak, method = "SJ")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Partlett-Riley prediction and confidence intervals
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: bias corrected Sidik-Jonkman
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3287 [-0.9835, 0.3261]
#>  d.f.: 8
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5625, -0.0950]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared:  85.5%
#> 
pimeta::pima(sbp$y, sbp$sigmak, method = "KR")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Partlett-Riley prediction and confidence intervals
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: Kenward-Roger
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3287 [-1.0280, 0.3706]
#>  d.f.: 5.9508
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5815, -0.0759]
#>  d.f.: 6.9508
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared:  85.5%
#> 
pimeta::pima(sbp$y, sbp$sigmak, method = "APX")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Partlett-Riley prediction and confidence intervals
#>  Heterogeneity variance: REML
#>  Variance for average treatment effect: approximate
#> 
#> No. of studies: 10
#> 
#> Average treatment effect [95% prediction interval]:
#>  -0.3287 [-0.9843, 0.3268]
#>  d.f.: 8
#> 
#> Average treatment effect [95% confidence interval]:
#>  -0.3287 [-0.5646, -0.0929]
#>  d.f.: 9
#> 
#> Heterogeneity measure
#>  tau-squared: 0.0700
#>  I-squared:  85.5%
#> 
```
