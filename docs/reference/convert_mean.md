# Converting means and standard deviations

Converting the means and standard deviations in the experimental and
control groups to the effect size estimates (mean differences) and the
within studies standard errors vector

## Usage

``` r
convert_mean(n1, m1, s1, n2, m2, s2, pooled = FALSE)
```

## Arguments

- n1:

  a vector of the number of patients in the experimental group

- m1:

  a vector of the means in the experimental group

- s1:

  a vector of the standard deviations in the experimental group

- n2:

  a vector of the number of patients in the control group

- m2:

  a vector of the means in the control group

- s2:

  a vector of the standard deviations in the control group

- pooled:

  logical; if `TRUE`, the pooled variance is used. The default is
  `FALSE`.

## Value

A `data.frame` of study data.

- `y`: a numeric vector of the effect size estimates (mean differences).

- `se`: a numeric vector of the within studies standard errors.

## Details

The effect size is the mean difference \\(y = m1 - m2\\). The within
study standard error is \\(\sqrt{s1^2/n1 + s2^2/n2}\\)
(`pooled = FALSE`, default) or \\(s_p \sqrt{1/n1 + 1/n2}\\) with the
pooled standard deviation \\(s_p^2 = ((n1-1) s1^2 + (n2-1)
s2^2)/(n1+n2-2)\\) (`pooled = TRUE`).

## See also

[`convert_bin`](https://nshi-stat.github.io/pimeta/reference/convert_bin.md),
[`clbp`](https://nshi-stat.github.io/pimeta/reference/clbp.md),
[`nfld`](https://nshi-stat.github.io/pimeta/reference/nfld.md)

## Examples

``` r
data(clbp, package = "pimeta")
dat <- pimeta::convert_mean(clbp$n1, clbp$m1, clbp$s1,
                            clbp$n2, clbp$m2, clbp$s2)
head(dat)
#>      y       se
#> 1  2.0 3.144921
#> 2  3.0 3.262923
#> 3 -2.0 4.132601
#> 4  3.3 3.745650
#> 5  5.0 3.745650
#> 6  8.3 5.884007
pimeta::pima(dat$y, dat$se, method = "HTS")
#> 
#> Prediction & Confidence Intervals for Random-Effects Meta-Analysis
#> 
#> Higgins-Thompson-Spiegelhalter prediction and confidence intervals
#>  Heterogeneity variance: DerSimonian-Laird
#>  Variance for average treatment effect: approximate
#> 
#> No. of studies: 23
#> 
#> Average treatment effect [95% prediction interval]:
#>  -3.1700 [-25.7503, 19.4103]
#>  d.f.: 21
#> 
#> Average treatment effect [95% confidence interval]:
#>  -3.1700 [-8.1178, 1.7778]
#>  d.f.: 22
#> 
#> Heterogeneity measure
#>  tau-squared: 112.2028
#>  I-squared:  92.2%
#> 
```
