# Converting binary data

Converting binary outcome data to the effect size estimates and the
within studies standard errors vector

## Usage

``` r
convert_bin(m1, n1, m2, n2, type = c("logOR", "logRR", "RD"))
```

## Arguments

- m1:

  a vector of the number of successes in the experimental group

- n1:

  a vector of the number of patients in the experimental group

- m2:

  a vector of the number of successes in the control group

- n2:

  a vector of the number of patients in the control group

- type:

  the outcome measure for binary outcome data (default = "logOR").

  - `logOR`: logarithmic odds ratio, which is defined by \\(=\log
    \frac{(m1+0.5)(n2-m2+0.5)}{(n1-m1+0.5)(m2+0.5)}\\).

  - `logRR`: logarithmic relative risk, which is defined by \\(=\log
    \frac{(m1+0.5)(n2+0.5)}{(n1+0.5)(m2+0.5)}\\).

  - `RD`: risk difference, which is defined by
    \\(=\frac{m1}{n1}-\frac{m2}{n2}\\).

## Value

A `data.frame` of study data.

- `y`: a numeric vector of the effect size estimates.

- `se`: a numeric vector of the within studies standard errors.

## Details

This function implements methods for logarithmic odds ratio, logarithmic
relative risk, and risk difference described in Hartung & Knapp (2001).

## References

Hartung, J., and Knapp, G. (2001). A refined method for the
meta-analysis of controlled clinical trials with binary outcome. *Stat
Med.* **20**(24): 3875-3889. <https://doi.org/10.1002/sim.1009>

## See also

[`convert_mean`](https://nshi-stat.github.io/pimeta/reference/convert_mean.md)

## Examples

``` r
m1 <- c(15,12,29,42,14,44,14,29,10,17,38,19,21)
n1 <- c(16,16,34,56,22,54,17,58,14,26,44,29,38)
m2 <- c( 9, 1,18,31, 6,17, 7,23, 3, 6,12,22,19)
n2 <- c(16,16,34,56,22,55,15,58,15,27,45,30,38)
dat <- pimeta::convert_bin(m1, n1, m2, n2, type = "logOR")
print(dat)
#>             y        se
#> 1   2.0989861 0.9847737
#> 2   3.3570262 1.0165653
#> 3   1.5652318 0.5747840
#> 4   0.8640463 0.4042977
#> 5   1.4656407 0.6332968
#> 6   2.2325713 0.4481371
#> 7   1.5465488 0.7782417
#> 8   0.4125323 0.3721812
#> 9   2.1202635 0.8265438
#> 10  1.8071598 0.6022988
#> 11  2.7646729 0.5382109
#> 12 -0.3544099 0.5555283
#> 13  0.2058521 0.4541130
pimeta::pima(dat$y, dat$se)
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> A parametric bootstrap prediction and confidence intervals
#>  Heterogeneity variance: DerSimonian-Laird
#>  Variance for average treatment effect: Hartung (Hartung-Knapp)
#> 
#> No. of studies: 13
#> 
#> Average treatment effect [95% prediction interval]:
#>  1.4209 [-0.6548, 3.5355]
#>  d.f.: 12
#> 
#> Average treatment effect [95% confidence interval]:
#>  1.4209 [0.7888, 2.0543]
#>  d.f.: 12
#> 
#> Heterogeneity measure
#>  tau-squared: 0.7176
#>  I-squared:  69.9%
#> 
```
