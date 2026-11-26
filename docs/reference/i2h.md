# \\(I^2\\) heterogeneity measure

Returns the estimator for (Higgins & Thompson, 2002).

## Usage

``` r
i2h(se, tau2h)
```

## Arguments

- se:

  the within studies standard errors vector

- tau2h:

  the estimate of \\(\tau^2\\)

## Value

- `i2h`: the estimate for \\(I^2\\).

## References

Higgins, J. P. T., and Thompson, S. G. (2002). Quantifying heterogeneity
in a meta-analysis. *Stat Med.* **21**(11): 1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## Examples

``` r
data(sbp, package = "pimeta")
tau2h <- pimeta::tau2h(sbp$y, sbp$sigmak)
pimeta::i2h(sbp$sigmak, tau2h$tau2h)
#> [1] 70.47668
```
