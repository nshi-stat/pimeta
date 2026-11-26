# Nonalcoholic fatty liver disease data

Data from the systematic review of Koutoukidis et al. (2019) of
randomised controlled trials of weight loss interventions in adults with
nonalcoholic fatty liver disease. Each row is a comparison between a
more-intensive weight loss intervention (experimental group) and no or a
lower-intensity weight loss intervention (control group); some studies
contribute more than one comparison. The outcome is the change in body
weight (kg) from baseline; negative values indicate weight loss. The
data can be converted to mean differences and their standard errors by
[`convert_mean`](https://nshi-stat.github.io/pimeta/reference/convert_mean.md).

## Usage

``` r
data(nfld)
```

## Format

A data frame with 25 rows and 8 variables

## Details

- `ID`: Comparison ID

- `Source`: First author name and year of publication

- `m1`: Estimated mean in the experimental group

- `s1`: Standard deviation in the experimental group

- `n1`: Number of observations in the experimental group

- `m2`: Estimated mean in the control group

- `s2`: Standard deviation in the control group

- `n2`: Number of observations in the control group

## References

Koutoukidis, D. A., Astbury, N. M., Tudor, K. E., Morris, E., Henry, J.
A., Noreik, M., Jebb, S. A., and Aveyard, P. (2019). Association of
weight loss interventions with changes in biomarkers of nonalcoholic
fatty liver disease: a systematic review and meta-analysis. *JAMA Intern
Med.* **179**(9): 1262-1271.
<https://doi.org/10.1001/jamainternmed.2019.2248>
