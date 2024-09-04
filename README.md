
<!-- README.md is generated from README.Rmd. Please edit that file -->

Selective reporting occurs when statistically significant, affirmative
results are more likely to be reported (and therefore more likely to be
available for meta-analysis) compared to null, non-affirmative results.
Selective reporting is a major concern for research synthesis because it
distorts the evidence base available for meta-analysis. Failure to
account for selective reporting can inflate effect size estimates from
meta-analysis and bias estimates of heterogeneity, making it difficult
to draw accurate conclusions from a synthesis.

There are many tools available already to investigate and correct for
selective outcome reporting. Widely used methods include: graphical
diagnostics like funnel plots; tests and adjustments for funnel plot
asymmetry like trim-and-fill, Egger’s regression, PET/PEESE, selection
models, and p-value diagnostics. However, very few methods available for
investigating selective reporting can accommodate dependent effect
sizes. Such limitation poses a problem for meta-analyses in education,
psychology and other social sciences, where dependent effects are a very
common feature of meta-analytic data.

Dependent effect sizes occur when primary studies report multiple
measures of the outcomes or repeated measures of the outcome. Failing to
account for dependency can result in misleading conclusions too narrow
confidence intervals and hypothesis tests that have inflated type one
error rates.

X (2024) developed and examined methods for investigating and accounting
for selective reporting in meta-analysis that account for dependent
effect sizes. Particularly, X (2024) combined step and beta selection
models with robust variance estimation and with cluster and fractional
weighted bootstrap. The results showed…

Our metaselection package provides a set of functions that implements
the methods discussed in X (2024).

## Installation

You can install the development version of the package from GitHub with:

``` r
remotes::install_github("jepusto/metaselection")
```

## Example

The following example uses `metadat::dat.lehmann` data from a
meta-analysis by Lehmann meta-analysis which examined the effects of
color red on attractiveness judgments. In the code below, we input the
`lehmann_dat` to the `selection_model()` function for our package to run
step function model with robust variance estimates.

``` r
library(metadat)
```

    ## Warning: package 'metadat' was built under R version 4.4.1

``` r
library(metaselection)

lehmann_dat <- dat.lehmann2018
lehmann_dat$sei <- sqrt(lehmann_dat$vi)

step_results <- selection_model(
  data = lehmann_dat, 
  yi = yi,
  sei = sei,
  selection_type = "step",
  steps = .025
)

step_results$est
```

    ##   estimator param        Est        SE      CI_lo      CI_hi
    ## 1        ML  beta  0.1327994 0.1377834 -0.1372511  0.4028498
    ## 2        ML gamma -2.5117243 1.0215649 -4.5139548 -0.5094939
    ## 3        ML zeta1 -0.6006530 1.0535870 -2.6656457  1.4643396

## Related Work

We want to recognize other packages that provide functions to selection
modeling.

The `metafor` package now includes the `selmodel()` function which
allows users to fit selection models. However, the function and the set
of selection models that it can fit can only be applied to meta-analytic
data assuming that the effects are independent.

## Acknowledgements

## References
