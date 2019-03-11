
# The 'pimeta' package


## Version 1.1.2 (2019-03-11)

* Prediction interval: the `pima` function
    - Parallel computing for the parametric bootstrap method (see a Vignette file).
    - Forest plot (see a Vignette file).
    - Kenward-Roger's approach.
* Confidence interval: the `cima` function
    - A Wald-type t-distribution confidence interval.
      Variance estimator of the average effect: an approximate estimator.
      Heterogeneity variance: Dearsimonian-Laird estimator.
    - A Wald-type t-distribution confidence interval.
      Variance estimator of the average effect: an approximate, Hartung-Knapp, Sidik-Jonkman, Kenward-Roger estimators.
      Heterogeneity variance: REML estimator.
    - Profile likelihood confidence interval.
    - Profile likelihood confidence interval with a Bartlett type correction.
    - Forest plot.
* Heterogeneity variance estimators: the `tau2h` function
    - DerSimonian-Laird estimator.
    - Variance component type estimator.
    - Paule--Mandel estimator.
    - Hartung-Makambi estimator.
    - Hunter--Schmidt estimator.
    - Maximum likelihood estimator.
    - Restricted maximum likelihood estimator.
    - Approximate restricted maximum likelihood estimator.
    - Sidik--Jonkman estimator.
    - Sidik--Jonkman improved estimator.
    - Empirical Bayes estimator.
    - Bayes modal estimator.
    - ML and REML confidence intervals.
* Converting binary data: the `convert_bin` function
    - Converting binary data to logarithmic odds ratio (see a Vignette file).
    - Converting binary data to logarithmic relative risk.
    - Converting binary data to risk difference.
* The distribution of a positive linear combination of chiqaure random variables: the `pwchisq` function


## Version 1.1.1 (2018-09-15)

* Fixed documents.


## Version 1.1.0 (2018-09-14)

* Refined the package structure.
* New function `pima` is available (see a Vignette file).


## Version 1.0.1 (2018-05-11)

* Updated citation informations.


## Version 1.0.0 (2018-04-05)

* First release.

