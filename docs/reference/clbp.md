# Chronic low back pain data

Data from the systematic review of Rubinstein et al. (2019) of
randomised controlled trials of spinal manipulative therapy for chronic
low back pain. Each row is a comparison between spinal manipulative
therapy (experimental group) and a recommended therapy (control group);
some studies contribute more than one comparison. The outcome is pain
intensity at one month on a 0-100 scale; lower values indicate less
pain. The data can be converted to mean differences and their standard
errors by
[`convert_mean`](https://nshi-stat.github.io/pimeta/reference/convert_mean.md).

## Usage

``` r
data(clbp)
```

## Format

A data frame with 23 rows and 8 variables

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

Rubinstein, S. M., de Zoete, A., van Middelkoop, M., Assendelft, W. J.
J., de Boer, M. R., and van Tulder, M. W. (2019). Benefits and harms of
spinal manipulative therapy for the treatment of chronic low back pain:
systematic review and meta-analysis of randomised controlled trials.
*BMJ.* **364**: l689.
[doi:10.1136/bmj.l689](https://doi.org/10.1136/bmj.l689)
