# Define prior penalty functions for selection model parameters

Creates a set of priors for use in estimating selection models. beta
parameters are assigned L-norm priors. log(tau) parameters are assigned
independent log-gamma priors. log(lambda) parameters are assigned
independent L-norm priors.

## Usage

``` r
define_priors(
  beta_mean = 0,
  beta_precision = 1/16,
  beta_L = 4,
  tau_mode = 0.2,
  tau_alpha = 1,
  lambda_mode = 0.5,
  lambda_precision = 1/54,
  lambda_L = 4
)
```

## Arguments

- beta_mean:

  numeric vector of prior means for beta (mean regression) parameters.

- beta_precision:

  numeric vector of prior precisions for beta (mean regression)
  parameters.

- beta_L:

  numeric vector of prior norms for beta (mean regression) parameters.

- tau_mode:

  numeric vector of prior modes for tau (heterogeneity SD) regression
  parameters.

- tau_alpha:

  numeric vector of prior precisions for tau (heterogeneity SD)
  regression parameters.

- lambda_mode:

  numeric vector of prior modes for lambda (selection) parameters.

- lambda_precision:

  numeric vector of prior precisions for lambda (selection) parameters.

- lambda_L:

  numeric vector of prior norms for lambda (mean regression) parameters.

## Value

An object of class `"selmodel_prior"` containing the following
components:

- `log_prior`:

  A function with arguments `beta`,`gamma`,`zeta` that returns the log
  of the prior density over these parameters.

- `score_prior`:

  A function with arguments `beta`,`gamma`,`zeta` that returns the
  vector of scores for the prior density over these parameters.

- `hessian_prior`:

  A function with arguments `beta`,`gamma`,`zeta` that returns the
  Hessian matrix of the prior density over these parameters.

## Examples

``` r
# set very informative priors on beta and lambda
strong_priors <- define_priors(
  beta_mean = 0.4, beta_precision = 40, 
  lambda_mode = 0.2, lambda_precision = 40
)

# set standard normal prior on beta
weak_priors <- define_priors(
  beta_mean = 0, beta_precision = 1/2, beta_L = 2
)
```
