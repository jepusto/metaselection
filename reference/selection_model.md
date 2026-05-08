# Estimate step or beta selection model

Estimate step or beta selection model, with standard errors and
confidence intervals based on either cluster-robust variance estimators
(i.e., sandwich estimators) or cluster-level bootstrapping to handle
dependent effect size estimates.

## Usage

``` r
selection_model(
  data,
  yi,
  vi,
  sei,
  pi,
  ai,
  cluster,
  selection_type = c("step", "beta"),
  steps = NULL,
  mean_mods = NULL,
  var_mods = NULL,
  sel_mods = NULL,
  sel_zero_mods = NULL,
  priors = define_priors(),
  subset = NULL,
  estimator = "CML",
  vcov_type = "robust",
  CI_type = "large-sample",
  conf_level = 0.95,
  theta = NULL,
  optimizer = NULL,
  optimizer_control = list(),
  use_jac = NULL,
  bootstrap = "none",
  R = 1999,
  retry_bootstrap = 0L,
  ...
)
```

## Arguments

- data:

  `data.frame` or `tibble` containing the meta-analytic data

- yi:

  vector of effect sizes estimates.

- vi:

  vector of sample variances. If `vi` is specified, then the `sei`
  argument must be omitted.

- sei:

  vector of sampling standard errors. If `sei` is specified, then the
  `vi` argument must be omitted.

- pi:

  optional vector of one-sided p-values. If not specified, p-values will
  be computed from `yi` and `sei`.

- ai:

  optional vector of analytic weights.

- cluster:

  vector indicating which observations belong to the same cluster.

- selection_type:

  character string specifying the type selection model to estimate, with
  possible options `"step"` or `"beta"`.

- steps:

  If `selection_type = "step"`, a numeric vector of one or more values
  specifying the thresholds (or steps) where the selection probability
  changes, with a default of `steps = .025`. If
  `selection_type = "beta"`, then a numeric vector of two values
  specifying the thresholds beyond which the selection function is
  truncated, with a default of `steps = c(.025, .975)`.

- mean_mods:

  optional model formula for moderators related to average effect size
  magnitude.

- var_mods:

  optional model formula for moderators related to effect size
  heterogeneity.

- sel_mods:

  optional model formula for moderators related to the probability of
  selection. Only relevant for `selection_type = "step"`.

- sel_zero_mods:

  optional model formula for moderators related to the probability of
  selection for p-values below the lowest threshold value of `steps`.
  Only relevant for `selection_type = "step"`.

- priors:

  a `selmodel_prior` object that defines priors (i.e., penalty terms)
  for model parameters, with a default of
  [`define_priors()`](define_priors.md). Set to `NULL` to obtain
  unpenalized estimates.

- subset:

  optional logical expression indicating a subset of observations to use
  for estimation.

- estimator:

  vector indicating whether to use the composite marginal likelihood
  estimator (option `"CML"`) or the augmented and reweighted Gaussian
  likelihood estimator (option `"ARGL"` or `"ARGL-full"`). If
  `selection_type = "beta"`, only the composite marginal likelihood
  estimator, `"CML"`, is available. For step function models, both
  estimators are available.

- vcov_type:

  character string specifying the type of variance-covariance matrix to
  calculate, with possible options `"robust"` for robust or
  cluster-robust standard errors, `"model-based"` for model-based
  standard errors, or `"none"`.

- CI_type:

  character string specifying the type of confidence interval to
  calculate, with possible options `"large-sample"` for large-sample
  normal interval (the default), `"percentile"` for a percentile
  interval, `"BCa"` for a bias-corrected-and-accelerated interval,
  `"bias-corrected"` for a bias-corrected percentile interval (without
  acceleration correction), `"normal"` for a standard normal interval,
  `"basic"` for a basic interval, `"student"` for a studentized
  interval, or `"none"`.

- conf_level:

  desired coverage level for confidence intervals, with the default
  value set to `.95`.

- theta:

  optional numeric vector of starting values to use in optimization
  routines.

- optimizer:

  character string indicating the optimizer to use. Ignored if
  `estimator = "ARGL"` or `"ARGL-full"`.

- optimizer_control:

  an optional list of control parameters to be used for optimization

- use_jac:

  logical indicating whether to use the Jacobian of the estimating
  equations for optimization. If `NULL` (the default), it will be reset
  to `FALSE` if `estimator = "CML"` or to `TRUE` if `estimator = "ARGL"`

- bootstrap:

  character string specifying the type of bootstrap to run, with
  possible options `"none"` (the default), `"exponential"` for the
  fractionally re-weighted cluster bootstrap, `"multinomial"` for a
  conventional clustered bootstrap, or , `"two-stage"` for a two-stage
  clustered bootstrap.

- R:

  number of bootstrap replications, with a default of `1999`.

- retry_bootstrap:

  number of times to re-draw a bootstrap sample in the event of
  non-convergence, with a default of `0`.

- ...:

  further arguments passed to
  [`simhelpers::bootstrap_CIs`](https://meghapsimatrix.github.io/simhelpers/reference/bootstrap_CIs.html).

## Value

An object of class `"selmodel"` containing the following components:

- `est`:

  A data frame with parameter estimates, standard errors, and confidence
  intervals. Note that the results do not include p-values so as to
  focus interpretation on the parameter estimates, rather than on the
  statistical significance of any given parameter.

- `vcov`:

  A matrix containing the estimated variance-covariance matrix of the
  parameter estiamtes

- `method`:

  Character string indicating the optimization method used to solve for
  parameter estimates.

- `info`:

  Further informaton about the optimization results.

- `ll`:

  Log likelihood of the model evaluated at the reported parameter
  estimates.

- `wpll`:

  Weighted partial log likelihood of the random effects model, with
  weights corresponding to inverse selection probabilities

- `n_clusters`:

  Number of independent clusters of effect sizes.

- `n_effects`:

  Number of effect size estimates in the data.

- `...`:

  Some additional elements containing information about the methods used
  to estimate the model.

## Examples

``` r
res_ML <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = 0.025,
  estimator = "CML",
  bootstrap = "none"
)

res_ML
#>    param    Est     SE  p_value   CI_lo CI_hi
#>     beta 0.2196 0.0525 2.85e-05 0.11677 0.322
#>     tau2 0.0394 0.0286       NA 0.00951 0.163
#>  lambda1 1.0350 0.5058 9.44e-01 0.39711 2.697
summary(res_ML)
#> Step Function Model 
#>  
#> Call: 
#> selection_model(data = self_control, yi = g, sei = se_g, cluster = studyid, 
#>     steps = 0.025, estimator = "CML", bootstrap = "none")
#> 
#> Number of clusters = 33; Number of effects = 158
#> 
#> Steps: 0.025 
#> Estimator: composite marginal likelihood 
#> Variance estimator: robust 
#> 
#> Log composite likelihood of selection model: -54.97404
#> Inverse selection weighted partial log likelihood: 87.11211 
#> 
#> Mean effect estimates:                                                
#>                                     Large Sample
#>  Coef. Estimate Std. Error  p-value Lower  Upper
#>   beta     0.22     0.0525 2.85e-05 0.117  0.322
#> 
#> Heterogeneity estimates:                                                 
#>                                      Large Sample
#>  Coef. Estimate Std. Error p-value   Lower  Upper
#>   tau2   0.0394     0.0286     --- 0.00951  0.163
#> 
#> Selection process estimates:
#>  Step: 0 < p <= 0.025; Studies: 18; Effects: 32                                                 
#>                                      Large Sample
#>    Coef. Estimate Std. Error p-value Lower  Upper
#>  lambda0        1        ---     ---   ---    ---
#> 
#>  Step: 0.025 < p <= 1; Studies: 24; Effects: 126                                                 
#>                                      Large Sample
#>    Coef. Estimate Std. Error p-value Lower  Upper
#>  lambda1     1.03      0.506   0.944 0.397    2.7

# configure progress bar
progressr::handlers(global = TRUE)
#> Error in globalCallingHandlers(condition = global_progression_handler): should not be called with handlers on the stack

res_hybrid <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = 0.025,
  estimator = "ARGL",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 19
)

res_hybrid
#>    param    Est     SE percentile_lower percentile_upper
#>     beta 0.2194 0.0483          0.07355            0.281
#>     tau2 0.0393 0.0296          0.00351            0.155
#>  lambda1 1.0332 0.1502          0.19581            2.674
summary(res_hybrid)
#> Step Function Model with Cluster Bootstrapping 
#>  
#> Call: 
#> selection_model(data = self_control, yi = g, sei = se_g, cluster = studyid, 
#>     steps = 0.025, estimator = "ARGL", CI_type = "percentile", 
#>     bootstrap = "multinomial", R = 19)
#> 
#> Number of clusters = 33; Number of effects = 158
#> 
#> Steps: 0.025 
#> Estimator: augmented and reweighted Gaussian likelihood 
#> Variance estimator: robust 
#> Bootstrap type: multinomial 
#> Number of bootstrap replications: 19 
#> 
#> Log composite likelihood of selection model: -54.9741
#> Inverse selection weighted partial log likelihood: 87.26052 
#> 
#> Mean effect estimates:                                               
#>                            Percentile Bootstrap
#>  Coef. Estimate Std. Error      Lower     Upper
#>   beta    0.219     0.0483     0.0735     0.281
#> 
#> Heterogeneity estimates:                                               
#>                            Percentile Bootstrap
#>  Coef. Estimate Std. Error      Lower     Upper
#>   tau2   0.0393     0.0296    0.00351     0.155
#> 
#> Selection process estimates:
#>  Step: 0 < p <= 0.025; Studies: 18; Effects: 32                                                 
#>                              Percentile Bootstrap
#>    Coef. Estimate Std. Error      Lower     Upper
#>  lambda0        1        ---        ---       ---
#> 
#>  Step: 0.025 < p <= 1; Studies: 24; Effects: 126                                                 
#>                              Percentile Bootstrap
#>    Coef. Estimate Std. Error      Lower     Upper
#>  lambda1     1.03       0.15      0.196      2.67
```
