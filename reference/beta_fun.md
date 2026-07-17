# Censor meta-analytic dataset based on the univariate beta-density model

A functional that takes model parameters and returns a function that can
be used to censor meta-analytic datasets according to the univariate
beta-density model.

## Usage

``` r
beta_fun(
  delta_1 = 1,
  delta_2 = 1,
  trunc_1 = 0.025,
  trunc_2 = 0.975,
  renormalize = TRUE
)
```

## Arguments

- delta_1:

  numeric value for the first parameter of the beta function

- delta_2:

  numeric value for the second parameter of the beta function

- trunc_1:

  numeric value between 0 and 1, below which p-values will be truncated.

- trunc_2:

  numeric value between 0 and 1, above which p-values will be truncated.

- renormalize:

  logical indicating whether to normalize the beta function to have a
  maximum value of 1, with a default value of `TRUE`.

## Value

A function that can be used to censor a meta-analytic dataset based on
the univariate beta-density model.
