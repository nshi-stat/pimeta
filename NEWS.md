
# The 'pimeta' package


## Version 1.1.4 (2026-11-26)

* Parametric bootstrap prediction interval (`method = "boot"`):
    - New argument `sampling`. The default `"table"` samples tau^2 from its
      confidence distribution through a quantile table (exact inversion at
      `nodes` probabilities and monotone cubic interpolation), which reduces
      the number of numerical inversions from `B` to about `nodes`.
      `sampling = "exact"` gives the previous per-sample inversion.
    - The R random number generator is no longer called inside the OpenMP
      parallel region; results with a given `seed` no longer depend on
      `parallel`.
    - `seed = NULL` no longer re-initializes the random number generator, so
      that `set.seed()` called before `pima()` gives reproducible results.
    - The `tol` argument of `pima()` was not passed to the bootstrap
      procedure; fixed.
    - A warning is issued when some bootstrap samples of tau^2 could not be
      generated.
    - The confidence distribution of tau^2 is computed on a standardised
      scale (`sigma / s`, `(y - weighted mean) / s`, `s^2 = mean(sigma^2)`).
      Previously `upper` and `tol` were absolute: when tau^2 exceeded 1000 on
      the original scale (e.g. mean differences in raw units) every sample
      was `Inf` and the interval was `NaN`, and on a very small scale `tol`
      exceeded tau^2 and the interval was far too wide. Results are now
      equivariant to the scale and location of `y`. The search is extended
      automatically (at least 1e8 times the mean of `sigma^2`).
    - Samples of tau^2 = `Inf` are handled as the limit tau^2 -> Inf (equal
      weights, predictive sample +-Inf) with a warning. Previously they gave
      `NaN` in the bootstrap samples (this happened for about 1 in 10^4
      samples for K = 3 with large heterogeneity).
    - When all `y` are identical, the sampling no longer loops (Q could be
      -1e-16 by rounding error); all samples of tau^2 are 0 and a warning
      is issued (the intervals have zero width).
    - Exact inversions outside the quantile table use the tolerance of the
      nodes (1e-10), and the search bracket is grown geometrically
      (faster).
* Heterogeneity variance estimators (the `tau2h` function):
    - `method = "VC"` always returned 0 (the sum of squared deviations was
      computed as the square of the sum, and `se` was used instead of
      `se^2`); fixed.
    - `method = "SJ"` and `"SJ2"` always returned (about) 0 because the
      residuals were not squared; fixed. `"SJ2"` is also affected by the
      fix of `"VC"` (its initial value).
    - `method = "HS"` used `Q - (K - 1)` instead of `Q - K` (Hunter and
      Schmidt, 2004; Viechtbauer, 2005); fixed.
    - The estimators used by `pima()` and `cima()` (DL, REML, ML) are not
      affected.
* Confidence interval (the `cima` function, `method = "boot"`): the `tol`
  argument was not passed to the bootstrap procedure; fixed. New arguments
  `sampling` and `nodes` (see `pima`).
* Fixed documentation errors (references, typos, and the variance of the
  average effect in the HTS method: `(sum w_i)^-1`).
* Added unit tests (testthat). The estimators are checked against
  independent implementations of the formulas and, if installed, against
  the metafor package.
* Forest plots (`plot.pima`, `plot.cima`): the deprecated
  `ggplot2::geom_errorbarh()` was replaced by `geom_linerange()`, and
  `size` by `linewidth`. ggplot2 (>= 3.4.0) is now required.
* Removed the unused `utils` from Imports.
* Examples no longer open a `cairo_pdf()` device.
* `inherits()` is used instead of comparing `class()` with a string.
* New function `convert_mean()`: mean differences and their standard
  errors from group means and standard deviations (unpooled or pooled
  variance).
* New datasets `clbp` (Rubinstein et al., 2019) and `nfld`
  (Koutoukidis et al., 2019).


## Version 1.1.3 (2019-09-17)

* Prediction & Confidence interval: the `pima` and `cima` functions
    - Back-transformartion for logarithmic scale outcomes (print & plot function; thanks to Dr. Morio Aihara).
* CITATION was updated.
* Vignette was updated.


## Version 1.1.2 (2019-03-11)

* Prediction interval: the `pima` function
    - Parallel computing for the parametric bootstrap method (see a Vignette file).
    - Forest plot (see a Vignette file).
    - Kenward-Roger's approach.
* Confidence interval: the `cima` function
    - A Wald-type t-distribution confidence interval.
      Variance estimator of the average effect: an approximate estimator.
      Heterogeneity variance: Dersimonian-Laird estimator.
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

