# Generate meta-analytic data

Generate meta-analytic correlated or correlated and hierarchical effects
data with options to simulate selective outcome reporting

## Usage

``` r
r_meta(
  mean_smd,
  tau,
  omega,
  m,
  cor_mu,
  cor_sd,
  n_ES_sim,
  censor_fun = NULL,
  m_multiplier = 1,
  id_start = 0L,
  paste_ids = TRUE,
  include_sel_prob = FALSE
)
```

## Arguments

- mean_smd:

  numeric value indicating the true mean effect size

- tau:

  numeric value characterizing between-study heterogeneity in effects

- omega:

  numeric value characterizing within-study heterogeneity in effects

- m:

  numeric value of studies in the simulated meta-analysis

- cor_mu:

  numeric value indicating the average correlation between outcomes

- cor_sd:

  numeric value indicating standard deviation of correlation between
  outcomes

- n_ES_sim:

  a function used to simulate the distribution of primary study sample
  sizes and the number of effect sizes per study

- censor_fun:

  a function used to censor effects; this package provides functionals
  [`step_fun()`](http://jepusto.github.io/metaselection/reference/step_fun.md)
  and
  [`beta_fun()`](http://jepusto.github.io/metaselection/reference/beta_fun.md)
  to censor effects based on step-function or beta-function models
  respectively. If `NULL` (the default)

- m_multiplier:

  numeric value indicating a multiplier for buffer for the number of
  studies

- id_start:

  integer indicating the starting value for study id

- paste_ids:

  logical with `TRUE` (the default) indicating that the study id and
  effect size id should be pasted together

- include_sel_prob:

  logical with `TRUE` indicating that the returned dataset should
  include a variable `selection_prob` reporting the true probability of
  selection given the observed p-value. Default of `FALSE` indicates
  that the `selection_prob` variable should be omitted.

## Value

A `data.frame` containing the simulated meta-analytic dataset.

## Examples

``` r

example_dat <- r_meta(
  mean_smd = 0,
  tau = .1, omega = .01,
  m = 50,
  cor_mu = .4, cor_sd = 0.001,
  censor_fun = step_fun(cut_vals = .025, weights = 0.4),
  n_ES_sim = n_ES_param(40, 3)
)
```
