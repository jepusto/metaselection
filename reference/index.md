# Package index

## Model-Fitting

Functions to fit and evaluate selection models with robust variance
estimation or cluster bootstrapping

- [`selection_model()`](http://jepusto.github.io/metaselection/reference/selection_model.md)
  : Estimate step or beta selection model

- [`selection_plot()`](http://jepusto.github.io/metaselection/reference/selection_plot.md)
  : Plot the selection weights implied by an estimated selection model.

- [`selection_wts()`](http://jepusto.github.io/metaselection/reference/selection_wts.md)
  : Calculate model-implied weights for specified p-values.

- [`print(`*`<selmodel>`*`)`](http://jepusto.github.io/metaselection/reference/print.selmodel.md)
  :

  Print results from a `selmodel` object

- [`summary(`*`<selmodel>`*`)`](http://jepusto.github.io/metaselection/reference/summary.selmodel.md)
  :

  Summarize results from a `selmodel` object

- [`define_priors()`](http://jepusto.github.io/metaselection/reference/define_priors.md)
  : Define prior penalty functions for selection model parameters

- [`p_area()`](http://jepusto.github.io/metaselection/reference/p_area.md)
  :

  Calculate area under the selection weight function from a `selmodel`
  object

## Example Datasets

Example datasets from published meta-analyses and one dataset containing
a distribution of primary study characteristics

- [`self_control`](http://jepusto.github.io/metaselection/reference/self_control.md)
  : Self-Control Training Meta-Analysis
- [`interleaved_learning`](http://jepusto.github.io/metaselection/reference/interleaved_learning.md)
  : Interleaved Learning Meta-Analysis
- [`practice_facilitation`](http://jepusto.github.io/metaselection/reference/practice_facilitation.md)
  : Practice Facilitation Meta-Analysis
- [`wwc_es`](http://jepusto.github.io/metaselection/reference/wwc_es.md)
  : What Works Clearinghouse sample size and effect size distribution
  data

## Data Simulation

Functions to simulate meta-analytic data, including data subject to
selective outcome reporting

- [`r_meta()`](http://jepusto.github.io/metaselection/reference/r_meta.md)
  : Generate meta-analytic data
- [`step_fun()`](http://jepusto.github.io/metaselection/reference/step_fun.md)
  : Censor meta-analytic dataset based on a univariate step-function
  model
- [`beta_fun()`](http://jepusto.github.io/metaselection/reference/beta_fun.md)
  : Censor meta-analytic dataset based on the univariate beta-density
  model
- [`step_count_fun()`](http://jepusto.github.io/metaselection/reference/step_count_fun.md)
  : Censor meta-analytic dataset based on a multivariate step-function
  model
- [`n_ES_empirical()`](http://jepusto.github.io/metaselection/reference/n_ES_empirical.md)
  : Simulate empirical distribution of sample size and number of effect
  sizes
- [`n_ES_param()`](http://jepusto.github.io/metaselection/reference/n_ES_param.md)
  : Simulate empirical distribution of sample size and number of effect
  sizes
