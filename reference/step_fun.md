# Censor meta-analytic dataset based on a univariate step-function model

A functional that takes in cut values and weights representing selection
probabilities for different intervals of p-values and returns a function
that can be used to censor meta-analytic datasets according to the
univariate step-function model.

## Usage

``` r
step_fun(cut_vals = 0.025, weights = 1, renormalize = TRUE)
```

## Arguments

- cut_vals:

  numeric vector of one or more values specifying the threshold (or
  step) where the selection probability changes.

- weights:

  numeric vector of one or more values specifying the selection
  probabilities for different intervals of p-values; the intervals are
  determined by the `cut_vals`.

- renormalize:

  logical indicating whether to normalize the step function to have a
  maximum value of 1, with a default value of `TRUE`.

## Value

A function that can be used to censor a meta-analytic dataset based on
the univariate step-function model.
