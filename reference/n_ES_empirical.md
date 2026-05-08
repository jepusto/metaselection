# Simulate empirical distribution of sample size and number of effect sizes

A functional that takes in a dataset with empirical distribution of
primary study sample sizes and number of effect sizes per primary study
and returns a function that generates random samples from the dataset

## Usage

``` r
n_ES_empirical(dat)
```

## Arguments

- dat:

  a `data.frame` or `tibble` containing primary study sample sizes and
  number of effect sizes per primary study

## Value

A function that generates random samples from the input dataset.

## Examples

``` r
study_features <- n_ES_empirical(wwc_es)
study_features(m=3)
#>       n n_ES
#> 173 605    4
#> 407 104    7
#> 460  80    2
study_features(m=7)
#>        n n_ES
#> 63  1458   22
#> 485   70    1
#> 559   48    3
#> 452   84    5
#> 201  467    9
#> 197  497    2
#> 517   61    1
```
