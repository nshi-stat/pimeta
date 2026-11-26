## Test environments

* local macOS Golden Gate 27.2.1 (aarch64-apple-darwin23), R 4.6.2
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes in this version

This is an update of an existing package (1.1.3 -> 1.1.4).

* Faster generation of tau^2 for the parametric bootstrap prediction
  interval (quantile table with monotone interpolation), with fixes for
  numerical problems at extreme scales.
* Bug fixes in several heterogeneity variance estimators (VC, SJ, SJ2,
  HS) in tau2h(); estimates from these methods differ from previous
  versions. The default methods used by pima() and cima() are not
  affected.
* New function convert_mean() and two example datasets (clbp, nfld).
* Added unit tests (testthat).
* Fixes for current R and ggplot2 (deprecated geom_errorbarh(),
  class() comparisons, examples using cairo_pdf()).

## Reverse dependencies

We checked 2 reverse dependencies, comparing R CMD check results across
CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

