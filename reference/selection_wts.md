# Calculate model-implied weights for specified p-values.

Calculates the selection weights implied by an estimated model for a
user-specified p-value or set of p-values.

## Usage

``` r
selection_wts(mod, pvals, ref_pval, ...)

# S3 method for class 'step.selmodel'
selection_wts(mod, pvals = NULL, ref_pval = NULL, bootstraps = TRUE, ...)

# S3 method for class 'beta.selmodel'
selection_wts(mod, pvals = NULL, ref_pval = NULL, bootstraps = TRUE, ...)
```

## Arguments

- mod:

  Fitted model of class `"selmodel"`.

- pvals:

  Numeric vector of p-values for which to calculate selection weights.

- ref_pval:

  Numeric value of a p-value at which to standardize the weights. If not
  `NULL`, then a p-value of `ref_pval` will have selection weight of 1
  and selection weights for all other p-values will be calculated
  relative to `ref_pval`.

- ...:

  further arguments passed to some methods.

- bootstraps:

  If `mod` includes bootstrap replications, then setting
  `bootstraps = TRUE` will return selection weights for each bootstrap
  replication, in addition to the selection weights implied by the model
  parameter estimates. Ignored if `mod` does not include bootstrap
  replications.

## Value

If `mod` does not include bootstrapped confidence intervals or if the
argument `bootstraps = FALSE`, then `selection_wts` will return a
`data.frame` containing the user-specified p-values and the selection
weights implied by the estimated model parameters.

If `mod` does include bootstrapped confidence intervals (i.e., when
`inherits(mod, "boot.selmodel")` is `TRUE`) and the argument
`bootstraps = TRUE`, then `selection_wts` will return a list with two
elements. The first element is a `data.frame` containing the
user-specified p-values and the selection weights implied by the
estimated model parameters. The second element is a `data.frame`
containing the user-specified p-values and the selection weights implied
by each bootstrap replicate of the model parameter estimates. The
`data.frame` includes an additional variable, `rep`, identifying the
bootstrap replicate.

## Examples

``` r
mod <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = c(0.025, .5),
  estimator = "CML",
  bootstrap = "none"
)

selection_wts(mod, pvals = seq(0, 1, 0.2))
#>     p        wt
#> 1 0.0 1.0000000
#> 2 0.2 0.8491666
#> 3 0.4 0.8491666
#> 4 0.6 0.3600605
#> 5 0.8 0.3600605
#> 6 1.0 0.3600605

mod_boot <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = c(0.025, .5),
  estimator = "CML",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 9
)

selection_wts(mod_boot, pvals = seq(0, 1, 0.2))
#> $wts
#>     p        wt
#> 1 0.0 1.0000000
#> 2 0.2 0.8491666
#> 3 0.4 0.8491666
#> 4 0.6 0.3600605
#> 5 0.8 0.3600605
#> 6 1.0 0.3600605
#> 
#> $boot_wts
#>    rep   p         wt
#> 1    1 0.0 1.00000000
#> 2    1 0.2 1.54116480
#> 3    1 0.4 1.54116480
#> 4    1 0.6 0.62523786
#> 5    1 0.8 0.62523786
#> 6    1 1.0 0.62523786
#> 7    2 0.0 1.00000000
#> 8    2 0.2 0.32455384
#> 9    2 0.4 0.32455384
#> 10   2 0.6 0.08847731
#> 11   2 0.8 0.08847731
#> 12   2 1.0 0.08847731
#> 13   3 0.0 1.00000000
#> 14   3 0.2 1.42939804
#> 15   3 0.4 1.42939804
#> 16   3 0.6 0.68465589
#> 17   3 0.8 0.68465589
#> 18   3 1.0 0.68465589
#> 19   4 0.0 1.00000000
#> 20   4 0.2 0.74693333
#> 21   4 0.4 0.74693333
#> 22   4 0.6 0.30575238
#> 23   4 0.8 0.30575238
#> 24   4 1.0 0.30575238
#> 25   5 0.0 1.00000000
#> 26   5 0.2 0.51347637
#> 27   5 0.4 0.51347637
#> 28   5 0.6 0.38449140
#> 29   5 0.8 0.38449140
#> 30   5 1.0 0.38449140
#> 31   6 0.0 1.00000000
#> 32   6 0.2 1.17760928
#> 33   6 0.4 1.17760928
#> 34   6 0.6 0.54865522
#> 35   6 0.8 0.54865522
#> 36   6 1.0 0.54865522
#> 37   7 0.0 1.00000000
#> 38   7 0.2 0.51054335
#> 39   7 0.4 0.51054335
#> 40   7 0.6 0.16658750
#> 41   7 0.8 0.16658750
#> 42   7 1.0 0.16658750
#> 43   8 0.0 1.00000000
#> 44   8 0.2 0.89145960
#> 45   8 0.4 0.89145960
#> 46   8 0.6 0.35845932
#> 47   8 0.8 0.35845932
#> 48   8 1.0 0.35845932
#> 49   9 0.0 1.00000000
#> 50   9 0.2 0.87269252
#> 51   9 0.4 0.87269252
#> 52   9 0.6 0.20720321
#> 53   9 0.8 0.20720321
#> 54   9 1.0 0.20720321
#> 

```
