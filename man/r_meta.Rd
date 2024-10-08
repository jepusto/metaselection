% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-generation.R
\name{r_meta}
\alias{r_meta}
\title{Generate meta-analytic data}
\usage{
r_meta(
  mean_smd,
  tau,
  omega,
  m,
  cor_mu,
  cor_sd,
  censor_fun,
  n_ES_sim,
  m_multiplier = 2,
  id_start = 0L,
  paste_ids = TRUE,
  include_sel_prob = FALSE
)
}
\arguments{
\item{mean_smd}{numeric value indicating the true mean effect size}

\item{tau}{numeric value characterizing between-study heterogeneity in
effects}

\item{omega}{numeric value characterizing within-study heterogeneity in
effects}

\item{m}{numeric value of studies in the simulated meta-analysis}

\item{cor_mu}{numeric value indicating the average correlation between
outcomes}

\item{cor_sd}{numeric value indicating standard deviation of correlation
between outcomes}

\item{censor_fun}{a function used to censor effects; this package provides
functionals \code{step_fun()} and \code{beta_fun()} to censor effects based on
step-function or beta-function models respectively.}

\item{n_ES_sim}{a function used to simulate the distribution of primary study
sample sizes and the number of effect sizes per study}

\item{m_multiplier}{numeric value indicating a multiplier for buffer for the
number of studies}

\item{id_start}{integer indicating the starting value for study id}

\item{paste_ids}{logical with \code{TRUE} (the default) indicating that the
study id and effect size id should be pasted together}

\item{include_sel_prob}{logical with \code{TRUE} indicating that the returned
dataset should include a variable \code{selection_prob} reporting the true
probability of selection given the observed p-value. Default of
\code{FALSE} indicates that the \code{selection_prob} variable should be
omitted.}
}
\value{
A \code{data.frame} containing the simulated meta-analytic dataset.
}
\description{
Generate meta-analytic correlated or correlated and hierarchical
effects data with options to simulate selective outcome reporting
}
\examples{

example_dat <- r_meta(
  mean_smd = 0,
  tau = .1, omega = .01,
  m = 50,
  cor_mu = .4, cor_sd = 0.001,
  censor_fun = step_fun(cut_vals = .025, weights = 0.4),
  n_ES_sim = n_ES_param(40, 3)
)

}
