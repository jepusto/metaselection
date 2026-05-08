# Simulate empirical distribution of sample size and number of effect sizes

A functional that takes in average sample size per primary study,
average number of effect sizes per study, and the minimum sample size
per study and returns a function to generate random samples of primary
study sample sizes and numbers of effect sizes per study.

## Usage

``` r
n_ES_param(mean_N, mean_ES, min_N = 20L)
```

## Arguments

- mean_N:

  numeric value specifying the average sample size per primary study

- mean_ES:

  numeric value specifying the average number of effect sizes per
  primary study

- min_N:

  numeric value specifying the minimum sample size per study

## Value

A function that generates a `data.frame` with randomly generated sample
size per primary study and number of effect sizes per study.

## Examples

``` r
study_features <- n_ES_param(mean_N = 40, mean_ES = 3, min_N = 10)
study_features(m = 3)
#>    n n_ES
#> 1 38    2
#> 2 38    4
#> 3 38    4
study_features(m = 5)
#>    n n_ES
#> 1 35    4
#> 2 31    6
#> 3 42    2
#> 4 35    2
#> 5 39    4
```
