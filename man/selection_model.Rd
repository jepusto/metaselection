% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_model.R
\name{selection_model}
\alias{selection_model}
\title{Estimate step or beta selection model}
\usage{
selection_model(
  data,
  yi,
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
  subset = NULL,
  estimator = "ML",
  vcov_type = "robust",
  CI_type = "large-sample",
  conf_level = 0.95,
  theta = NULL,
  optimizer = NULL,
  optimizer_control = list(),
  use_jac = TRUE,
  bootstrap = "none",
  R = 1999,
  ...
)
}
\arguments{
\item{data}{\code{data.frame} or \code{tibble} containing the meta-analytic
data}

\item{yi}{vector of effect sizes estimates.}

\item{sei}{vector of sampling standard errors.}

\item{pi}{optional vector of one-sided p-values. If not specified, p-values
will be computed from \code{yi} and \code{sei}.}

\item{ai}{optional vector of analytic weights.}

\item{cluster}{vector indicating which observations belong to the same
cluster.}

\item{selection_type}{character string specifying the type selection model to
estimate, with possible options \code{"step"} or \code{"beta"}.}

\item{steps}{If \code{selection_type = "step"}, a numeric vector of one or
more values specifying the thresholds (or steps) where the selection
probability changes, with a default of \code{steps = .025}. If
\code{selection_type = "beta"}, then a numeric vector of two values
specifying the thresholds beyond which the selection function is truncated,
with a default of \code{steps = c(.025, .975)}.}

\item{mean_mods}{optional model formula for moderators related to average
effect size magnitude.}

\item{var_mods}{optional model formula for moderators related to effect size
heterogeneity.}

\item{sel_mods}{optional model formula for moderators related to the
probability of selection. Only relevant for \code{selection_type = "step"}.}

\item{sel_zero_mods}{optional model formula for moderators related to the
probability of selection for p-values below the lowest threshold value of
\code{steps}. Only relevant for \code{selection_type = "step"}.}

\item{subset}{optional logical expression indicating a subset of observations
to use for estimation.}

\item{estimator}{vector indicating whether to use the maximum likelihood or
the hybrid estimator, with possible options \code{"ML"}, \code{"hybrid"},
and \code{"hybrid-full"}. If \code{selection_type = "beta"}, only the
maximum likelihood estimator, \code{"ML"}, is available. For step function
models, both maximum likelihood and hybrid estimators are available.}

\item{vcov_type}{character string specifying the type of variance-covariance
matrix to calculate, with possible options \code{"robust"} for robust or
cluster-robust standard errors, \code{"model-based"} for model-based
standard errors, or \code{"none"}.}

\item{CI_type}{character string specifying the type of confidence interval to
calculate, with possible options \code{"large-sample"} for large-sample
normal (the default), \code{"percentile"} for percentile bootstrap,
\code{"basic"} for basic bootstrap, \code{"student"} for studentized
bootstrap, or \code{"none"}.}

\item{conf_level}{desired coverage level for confidence intervals, with the
default value set to \code{.95}}

\item{theta}{optional numeric vector of starting values to use in
optimization routines.}

\item{optimizer}{character string indicating the optimizer to use. Ignored if
\code{estimator = "hybrid"} or \code{"hybrid-full"}.}

\item{optimizer_control}{an optional list of control parameters to be used
for optimization}

\item{use_jac}{logical with \code{TRUE} (the default) indicating to use the
Jacobian of the estimating equations for optimization.}

\item{bootstrap}{character string specifying the type of bootstrap to run,
with possible options \code{"none"} (the default), \code{"exponential"} for
the fractionally re-weighted cluster bootstrap, or \code{"multinomial"} for
a conventional clustered bootstrap.}

\item{R}{number of bootstrap replications, with a default of \code{1999}.}

\item{...}{further arguments passed to \code{simhelpers::bootstrap_CIs}.}
}
\value{
A numeric vector.
}
\description{
Estimate step or beta selection model
}
\examples{
res_ML <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = 0.025,
  estimator = "ML",
  bootstrap = "none"
)

res_ML
summary(res_ML)

# configure progress bar
progressr::handlers(global = TRUE)

res_hybrid <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = 0.025,
  estimator = "hybrid",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 19
)

res_hybrid
summary(res_hybrid)

}
