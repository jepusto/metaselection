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
#> 2    1 0.2 0.42079907
#> 3    1 0.4 0.42079907
#> 4    1 0.6 0.57510991
#> 5    1 0.8 0.57510991
#> 6    1 1.0 0.57510991
#> 7    2 0.0 1.00000000
#> 8    2 0.2 0.74076383
#> 9    2 0.4 0.74076383
#> 10   2 0.6 0.42729261
#> 11   2 0.8 0.42729261
#> 12   2 1.0 0.42729261
#> 13   3 0.0 1.00000000
#> 14   3 0.2 1.09165960
#> 15   3 0.4 1.09165960
#> 16   3 0.6 0.67490877
#> 17   3 0.8 0.67490877
#> 18   3 1.0 0.67490877
#> 19   4 0.0 1.00000000
#> 20   4 0.2 0.98564070
#> 21   4 0.4 0.98564070
#> 22   4 0.6 0.37402561
#> 23   4 0.8 0.37402561
#> 24   4 1.0 0.37402561
#> 25   5 0.0 1.00000000
#> 26   5 0.2 0.26783717
#> 27   5 0.4 0.26783717
#> 28   5 0.6 0.06796321
#> 29   5 0.8 0.06796321
#> 30   5 1.0 0.06796321
#> 31   6 0.0 1.00000000
#> 32   6 0.2 0.55984146
#> 33   6 0.4 0.55984146
#> 34   6 0.6 0.18824126
#> 35   6 0.8 0.18824126
#> 36   6 1.0 0.18824126
#> 37   7 0.0 1.00000000
#> 38   7 0.2 0.67724513
#> 39   7 0.4 0.67724513
#> 40   7 0.6 0.43466100
#> 41   7 0.8 0.43466100
#> 42   7 1.0 0.43466100
#> 43   8 0.0 1.00000000
#> 44   8 0.2 0.90167691
#> 45   8 0.4 0.90167691
#> 46   8 0.6 0.51721245
#> 47   8 0.8 0.51721245
#> 48   8 1.0 0.51721245
#> 49   9 0.0 1.00000000
#> 50   9 0.2 0.24975183
#> 51   9 0.4 0.24975183
#> 52   9 0.6 0.18292967
#> 53   9 0.8 0.18292967
#> 54   9 1.0 0.18292967
#> 

```
