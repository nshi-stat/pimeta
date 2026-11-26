# Plot Results

A function for plotting of \`pima\` objects.

## Usage

``` r
# S3 method for class 'pima'
plot(
  x,
  y = NULL,
  title = "Forest plot",
  base_size = 16,
  base_family = "",
  digits = 3,
  studylabel = NULL,
  ntick = NULL,
  trans = c("identity", "exp"),
  ...
)
```

## Arguments

- x:

  \`pima\` object to plot

- y:

  is not used

- title:

  graph title

- base_size:

  base font size

- base_family:

  base font family

- digits:

  a value for digits specifies the minimum number of significant digits
  to be printed in values.

- studylabel:

  labels for each study

- ntick:

  the number of x-axis ticks

- trans:

  transformation for logarithmic scale outcomes (`"identity"`
  \[default\] or `"exp"`).

- ...:

  further arguments passed to or from other methods.

## Examples

``` r
# \donttest{
data(sbp, package = "pimeta")
piex <- pimeta::pima(sbp$y, sbp$sigmak, method = "HTS")
plot(piex, digits = 2, base_size = 10, studylabel = sbp$label)

# }
```
