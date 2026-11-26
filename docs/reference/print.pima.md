# Print Results

`print` prints its argument and returns it invisibly (via
`invisible(x)`).

## Usage

``` r
# S3 method for class 'pima'
print(x, digits = 4, trans = c("identity", "exp"), ...)
```

## Arguments

- x:

  print to display

- digits:

  a value for digits specifies the minimum number of significant digits
  to be printed in values.

- trans:

  transformation for logarithmic scale outcomes (`"identity"`
  \[default\] or `"exp"`).

- ...:

  further arguments passed to or from other methods.
