# Calculate area under the selection weight function from a `selmodel` object

Summarize the strength of selection by calculating the area under the
selection weight function a `selmodel` object (excluding the area with
weight fixed at 1). If the object has bootstrap replications, then a
confidence interval will also be calculated.

## Usage

``` r
p_area(object, CI_type = NULL, conf_level = NULL, warn = TRUE)
```

## Arguments

- object:

  Fitted model of class `"selmodel"`.

- CI_type:

  character string specifying the type of confidence interval to
  calculate, with options as in `"selection_model"`. If `NULL` (the
  default), it will be inherited from `object`.

- conf_level:

  desired coverage level for confidence intervals. If `NULL` (the
  default), it will be inherited from `object`, which has a default
  value of `.95`.

- warn:

  logical controlling whether warnings are displayed, with a default of
  `TRUE`.

## Examples

``` r

beta_noboot <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  selection_type = "beta",
  steps = c(0.025,0.5)
)

p_area(beta_noboot)
#>    param       Est
#> 1 p-area 0.6420617

step_boot <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  selection_type = "step",
  steps = c(0.025,0.50),
  estimator = "ARGL",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 6
)

p_area(step_boot)
#>    param       Est       SE bootstraps percentile_lower percentile_upper
#> 1 p-area 0.5989296 5.284236          6       0.09451833         13.41982
```
