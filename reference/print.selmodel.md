# Print results from a `selmodel` object

Print relevant results from a fitted `selmodel` object.

## Usage

``` r
# S3 method for class 'selmodel'
print(x, transf_gamma = TRUE, transf_zeta = TRUE, digits = 3, ...)
```

## Arguments

- x:

  Fitted model of class `"selmodel"`.

- transf_gamma:

  logical with `TRUE` (the default) indicating that the heterogeneity
  parameter estimates (called gamma) should be transformed by
  exponentiating.

- transf_zeta:

  logical with `TRUE` (the default) indicating that the selection
  parameter estimates (called zeta) should be transformed by
  exponentiating.

- digits:

  Minimum number of significant digits to be used, with a default of 3.

- ...:

  further arguments passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

## Examples

``` r
res_ML <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = 0.025,
  estimator = "CML",
  bootstrap = "none"
)

print(res_ML)
#>    param    Est     SE  p_value   CI_lo CI_hi
#>     beta 0.2196 0.0525 2.85e-05 0.11677 0.322
#>     tau2 0.0394 0.0286       NA 0.00951 0.163
#>  lambda1 1.0350 0.5058 9.44e-01 0.39711 2.697
print(res_ML, transf_gamma = FALSE, transf_zeta = FALSE)
#>  param     Est     SE  p_value  CI_lo  CI_hi
#>   beta  0.2196 0.0525 2.85e-05  0.117  0.322
#>  gamma -3.2338 0.7252       NA -4.655 -1.812
#>  zeta1  0.0344 0.4887 9.44e-01 -0.924  0.992
 
```
