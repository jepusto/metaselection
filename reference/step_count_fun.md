# Censor meta-analytic dataset based on a multivariate step-function model

A functional that takes in a single step value, a weight representing
the selection probability for the upper interval of p-values, and a
dependence parameter, and returns a function that can be used to censor
meta-analytic datasets according to a multivariate step-function model.

## Usage

``` r
step_count_fun(cut_val = 0.025, weight = 1, psi = 0, renormalize = TRUE)
```

## Arguments

- cut_val:

  numeric value specifying the specifying the thresholds (or steps)
  where the selection probability changes.

- weight:

  numeric value specifying the selection probability for the upper
  interval of p-values, i.e., for p-values larger than `cut_val` for an
  effect reported in a study with no p-values smaller than `cut_val`.

- psi:

  numeric value controlling the degree of dependence in selection
  probabilities.

- renormalize:

  logical indicating whether to normalize the step function to have a
  maximum value of 1, with a default value of `TRUE`.

## Value

A function that can be used to censor a meta-analytic dataset based on a
multivariate step-function model.
