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

  fitted model of class `"selmodel"`.

- pvals:

  numeric vector of p-values for which to calculate selection weights.

- ref_pval:

  numeric value of a p-value at which to standardize the weights. If not
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
  estimator = "ARGL"
)

selection_wts(mod, pvals = seq(0, 1, 0.2))
#>     p        wt
#> 1 0.0 1.0000000
#> 2 0.2 0.8490364
#> 3 0.4 0.8490364
#> 4 0.6 0.3613281
#> 5 0.8 0.3613281
#> 6 1.0 0.3613281

mod_boot <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = c(0.025, .5),
  estimator = "ARGL",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 9
)

selection_wts(mod_boot, pvals = seq(0, 1, 0.2))
#> $wts
#>     p        wt
#> 1 0.0 1.0000000
#> 2 0.2 0.8490364
#> 3 0.4 0.8490364
#> 4 0.6 0.3613281
#> 5 0.8 0.3613281
#> 6 1.0 0.3613281
#> 
#> $boot_wts
#>    rep   p         wt
#> 1    1 0.0 1.00000000
#> 2    1 0.2 2.11960071
#> 3    1 0.4 2.11960071
#> 4    1 0.6 2.04984288
#> 5    1 0.8 2.04984288
#> 6    1 1.0 2.04984288
#> 7    2 0.0 1.00000000
#> 8    2 0.2 0.63116523
#> 9    2 0.4 0.63116523
#> 10   2 0.6 0.25100172
#> 11   2 0.8 0.25100172
#> 12   2 1.0 0.25100172
#> 13   3 0.0 1.00000000
#> 14   3 0.2 0.77718531
#> 15   3 0.4 0.77718531
#> 16   3 0.6 0.26483471
#> 17   3 0.8 0.26483471
#> 18   3 1.0 0.26483471
#> 19   4 0.0 1.00000000
#> 20   4 0.2 0.25538886
#> 21   4 0.4 0.25538886
#> 22   4 0.6 0.05548891
#> 23   4 0.8 0.05548891
#> 24   4 1.0 0.05548891
#> 25   5 0.0 1.00000000
#> 26   5 0.2 0.94684842
#> 27   5 0.4 0.94684842
#> 28   5 0.6 0.23342662
#> 29   5 0.8 0.23342662
#> 30   5 1.0 0.23342662
#> 31   6 0.0 1.00000000
#> 32   6 0.2 0.44373011
#> 33   6 0.4 0.44373011
#> 34   6 0.6 0.34469000
#> 35   6 0.8 0.34469000
#> 36   6 1.0 0.34469000
#> 37   7 0.0 1.00000000
#> 38   7 0.2 0.35494918
#> 39   7 0.4 0.35494918
#> 40   7 0.6 0.22472733
#> 41   7 0.8 0.22472733
#> 42   7 1.0 0.22472733
#> 43   8 0.0 1.00000000
#> 44   8 0.2 1.74742269
#> 45   8 0.4 1.74742269
#> 46   8 0.6 0.82046800
#> 47   8 0.8 0.82046800
#> 48   8 1.0 0.82046800
#> 49   9 0.0 1.00000000
#> 50   9 0.2 4.71042395
#> 51   9 0.4 4.71042395
#> 52   9 0.6 3.46602561
#> 53   9 0.8 3.46602561
#> 54   9 1.0 3.46602561
#> 

```
