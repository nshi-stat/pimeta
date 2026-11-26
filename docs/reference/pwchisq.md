# The Distribution of a Positive Linear Combination of Chiqaure Random Variables

The cumulative distribution function for the distribution of a positive
linear combination of \\(\chi^2\\) random variables with the weights
(\\(\lambda_1, \ldots, \lambda_K\\)), degrees of freedom (\\(\nu_1,
\ldots, \nu_K\\)), and non-centrality parameters (\\(\delta_1, \ldots,
\delta_K\\)).

## Usage

``` r
pwchisq(
  x,
  lambda = 1,
  nu = 1,
  delta = 0,
  mode = 1,
  maxit1 = 1e+05,
  eps = 10^(-10)
)
```

## Arguments

- x:

  numeric; value of x \> 0 (\\(P\[X \le x\]\\)).

- lambda:

  numeric vector; weights (\\(\lambda_1, \ldots, \lambda_K\\)).

- nu:

  integer vector; degrees of freedom (\\(\nu_1, \ldots, \nu_K\\)).

- delta:

  numeric vector; non-centrality parameters (\\(\delta_1, \ldots,
  \delta_K\\)).

- mode:

  numeric; the mode of calculation (see Farabrother, 1984)

- maxit1:

  integer; the maximum number of iteration.

- eps:

  numeric; the desired level of accuracy.

## Value

- `prob`: the distribution function.

## References

Farebrother, R. W. (1984). Algorithm AS 204: the distribution of a
positive linear combination of \\(\chi^2\\) random variables. *J R Stat
Soc Ser C Appl Stat.* **33**(3): 332–339.
<https://rss.onlinelibrary.wiley.com/doi/10.2307/2347721>.

## Examples

``` r
# Table 1 of Farabrother (1984)
# Q6 (The taget values are 0.0061, 0.5913, and 0.9779)

pimeta::pwchisq( 20, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#> [1] 0.006117973
pimeta::pwchisq(100, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#> [1] 0.5913421
pimeta::pwchisq(200, lambda = c(7,3), nu = c(6,2), delta = c(6,2))
#> [1] 0.9779184
# [1] 0.006117973
# [1] 0.5913421
# [1] 0.9779184
```
