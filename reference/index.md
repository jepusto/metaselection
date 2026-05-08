# Package index

## Model-Fitting

Functions to fit and evaluate selection models with robust variance
estimation or cluster bootstrapping

- [`selection_model()`](selection_model.md) : Estimate step or beta
  selection model

- [`selection_plot()`](selection_plot.md) : Plot the selection weights
  implied by an estimated selection model.

- [`selection_wts()`](selection_wts.md) : Calculate model-implied
  weights for specified p-values.

- [`print(`*`<selmodel>`*`)`](print.selmodel.md) :

  Print results from a `selmodel` object

- [`summary(`*`<selmodel>`*`)`](summary.selmodel.md) :

  Summarize results from a `selmodel` object

- [`define_priors()`](define_priors.md) : Define prior penalty functions
  for selection model parameters

- [`p_area()`](p_area.md) :

  Calculate area under the selection weight function from a `selmodel`
  object

## Example Datasets

Example datasets from published meta-analyses and one dataset containing
a distribution of primary study characteristics

- [`self_control`](self_control.md) : Self-Control Training
  Meta-Analysis
- [`interleaved_learning`](interleaved_learning.md) : Interleaved
  Learning Meta-Analysis
- [`practice_facilitation`](practice_facilitation.md) : Practice
  Facilitation Meta-Analysis
- [`wwc_es`](wwc_es.md) : What Works Clearinghouse sample size and
  effect size distribution data

## Data Simulation

Functions to simulate meta-analytic data, including data subject to
selective outcome reporting

- [`r_meta()`](r_meta.md) : Generate meta-analytic data
- [`step_fun()`](step_fun.md) : Censor meta-analytic dataset based on a
  step-function model
- [`beta_fun()`](beta_fun.md) : Censor meta-analytic dataset based on
  the beta-density model
- [`n_ES_empirical()`](n_ES_empirical.md) : Simulate empirical
  distribution of sample size and number of effect sizes
- [`n_ES_param()`](n_ES_param.md) : Simulate empirical distribution of
  sample size and number of effect sizes
