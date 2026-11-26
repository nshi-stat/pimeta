# Calculating Heterogeneity Variance

Returns a heterogeneity variance estimate and its confidence interval.

## Usage

``` r
tau2h(
  y,
  se,
  maxiter = 100,
  method = c("DL", "VC", "PM", "HM", "HS", "ML", "REML", "AREML", "SJ", "SJ2", "EB",
    "BM"),
  methodci = c(NA, "ML", "REML"),
  alpha = 0.05
)
```

## Arguments

- y:

  the effect size estimates vector

- se:

  the within studies standard errors vector

- maxiter:

  the maximum number of iterations

- method:

  the calculation method for heterogeneity variance (default = "DL").

  - `DL`: DerSimonian–Laird estimator (DerSimonian & Laird, 1986).

  - `VC`: Variance component type estimator (Hedges, 1983).

  - `PM`: Paule–Mandel estimator (Paule & Mandel, 1982).

  - `HM`: Hartung–Makambi estimator (Hartung & Makambi, 2003).

  - `HS`: Hunter–Schmidt estimator (Hunter & Schmidt, 2004). This
    estimator has negative bias (Viechtbauer, 2005).

  - `ML`: Maximum likelihood (ML) estimator (e.g., DerSimonian & Laird,
    1986).

  - `REML`: Restricted maximum likelihood (REML) estimator (e.g.,
    DerSimonian & Laird, 1986).

  - `AREML`: Approximate restricted maximum likelihood estimator
    (Thompson & Sharp, 1999).

  - `SJ`: Sidik–Jonkman estimator (Sidik & Jonkman, 2005).

  - `SJ2`: Sidik–Jonkman improved estimator (Sidik & Jonkman, 2007).

  - `EB`: Empirical Bayes estimator (Morris, 1983).

  - `BM`: Bayes modal estimator (Chung, et al., 2013).

- methodci:

  the calculation method for a confidence interval of heterogeneity
  variance (default = NA).

  - `NA`: a confidence interval will not be calculated.

  - `ML`: Wald confidence interval with a ML estimator (Biggerstaff &
    Tweedie, 1997).

  - `REML`: Wald confidence interval with a REML estimator (Biggerstaff
    & Tweedie, 1997).

- alpha:

  the alpha level of the confidence interval

## Value

- `tau2h`: the estimate for \\(\tau^2\\).

- `lci`, `uci`: the lower and upper confidence limits
  \\(\hat{\tau}^2_l\\) and \\(\hat{\tau}^2_u\\).

## Details

Excellent reviews of heterogeneity variance estimation have been
published (Sidik & Jonkman, 2007; Veroniki, et al., 2016; Langan, et
al., 2019).

## References

Sidik, K., and Jonkman, J. N. (2007). A comparison of heterogeneity
variance estimators in combining results of studies. *Stat Med.*
**26**(9): 1964-1981. <https://doi.org/10.1002/sim.2688>.

Veroniki, A. A., Jackson, D., Viechtbauer, W., Bender, R., Bowden, J.,
Knapp, G., Kuss, O., Higgins, J. P. T., Langan, D., and Salanti, J.
(2016). Methods to estimate the between-study variance and its
uncertainty in meta-analysis. *Res Syn Meth.* **7**(1): 55-79.
<https://doi.org/10.1002/jrsm.1164>.

Langan, D., Higgins, J. P. T., Jackson, D., Bowden, J., Veroniki, A. A.,
Kontopantelis, E., Viechtbauer, W., and Simmonds, M. (2019). A
comparison of heterogeneity variance estimators in simulated
random-effects meta-analyses. *Res Synth Methods.* **10**(1): 83-98.
<https://doi.org/10.1002/jrsm.1316>.

DerSimonian, R., and Laird, N. (1986). Meta-analysis in clinical trials.
*Control Clin Trials.* **7**(3): 177-188.
<https://doi.org/10.1016/0197-2456(86)90046-2>.

Hedges, L. V. (1983). A random effects model for effect sizes. *Psychol
Bull.* **93**(2): 388-395. <https://doi.org/10.1037/0033-2909.93.2.388>.

Paule, R. C., and Mandel, K. H. (1982). Consensus values and weighting
factors. *J Res Natl Inst Stand Techno.* **87**(5): 377-385.
<https://doi.org/10.6028/jres.087.022>.

Hartung, J., and Makambi, K. H. (2003). Reducing the number of
unjustified significant results in meta-analysis. *Commun Stat Simul
Comput.* **32**(4): 1179-1190. <https://doi.org/10.1081/SAC-120023884>.

Hunter, J. E., and Schmidt, F. L. (2004). *Methods of Meta-Analysis:
Correcting Error and Bias in Research Findings. 2nd edition.* Sage
Publications, Inc.

Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance
estimators in the random-effects model. *J Educ Behav Stat.* **30**(3):
261-293. <https://doi.org/10.3102/10769986030003261>.

Thompson, S. G., and Sharp, S. J. (1999). Explaining heterogeneity in
meta-analysis: a comparison of methods. *Stat Med.* **18**(20):
2693-2708.
<https://doi.org/10.1002/(SICI)1097-0258(19991030)18:20%3C2693::AID-SIM235%3E3.0.CO;2-V>.

Sidik, K., and Jonkman, J. N. (2005). Simple heterogeneity variance
estimation for meta-analysis. *J R Stat Soc Ser C Appl Stat.* **54**(2):
367-384. <https://doi.org/10.1111/j.1467-9876.2005.00489.x>.

Morris, C. N. (1983). Parametric empirical Bayes inference: theory and
applications. *J Am Stat Assoc.* **78**(381): 47-55.
<https://doi.org/10.1080/01621459.1983.10477920>.

Chung, Y. L., Rabe-Hesketh, S., and Choi, I-H. (2013). Avoiding zero
between-study variance estimates in random-effects meta-analysis. *Stat
Med.* **32**(23): 4071-4089. <https://doi.org/10.1002/sim.5821>.

Biggerstaff, B. J., and Tweedie, R. L. (1997). Incorporating variability
in estimates of heterogeneity in the random effects model in
meta-analysis. *Stat Med.* **16**(7): 753-768.
<https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7%3C753::AID-SIM494%3E3.0.CO;2-G>.

## Examples

``` r
data(sbp, package = "pimeta")
pimeta::tau2h(sbp$y, sbp$sigmak)
#> 
#> Heterogeneity Variance for Random-Effects Meta-Analysis
#> 
#> Estimation method(s)
#>  Point estimation: Dersimonian-Laird estimator
#> 
#> No. of studies: 10
#> 
#> tau-squared: 0.028
#> 
```
