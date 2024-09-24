
<!-- badges: start -->

[![R-CMD-check](https://github.com/jepusto/metaselection/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jepusto/metaselection/actions/workflows/R-CMD-check.yaml)
[![Codecov
Status](https://codecov.io/gh/jepusto/metaselection/graph/badge.svg?token=8T7IUFT1QV)](https://codecov.io/gh/jepusto/metaselection)
<!-- [![CRAN Version](http://www.r-pkg.org/badges/version/metaselection)](https://CRAN.R-project.org/package=metaselection) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/grand-total/metaselection)](https://CRAN.R-project.org/package=metaselection) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/last-month/metaselection)](https://CRAN.R-project.org/package=metaselection) -->
<!-- badges: end -->

# metaselection

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
effect sizes. The results showed that selection models combined with
robust variance estimation led to lower bias in the estimate of the
overall effect size. Combining the selection models with cluster
bootstrapping led to close to nominal coverage rates.

Our metaselection package provides a set of functions that implements
these methods. The main function, `selection_model()`, fits step and
beta selection models with robust variance estimation and has options to
run cluster bootstrapping.

## Installation

You can install the development version of the package from GitHub with:

``` r
remotes::install_github("jepusto/metaselection")
```

## Example

The following example uses `metadat::dat.lehmann` data from a
meta-analysis by Lehmann meta-analysis which examined the effects of
color red on attractiveness judgments (White et al. 2022). In the code
below, we input the `lehmann_dat` to the `selection_model()` function
for our package to run step function model with cluster bootstrapping.
For further details, please see the vignette.

``` r
library(metaselection)

data("dat.lehmann2018", package = "metadat")
dat.lehmann2018$study <- dat.lehmann2018$Full_Citation
dat.lehmann2018$sei <- sqrt(dat.lehmann2018$vi)

set.seed(20240910)

mod_3PSM_boot <- selection_model(
  data = dat.lehmann2018, 
  yi = yi,
  sei = sei,
  cluster = study,
  selection_type = "step",
  steps = .025,
  CI_type = "percentile",
  bootstrap = "multinomial",
  R = 19
)

print(mod_3PSM_boot, transf_gamma = TRUE, transf_zeta = TRUE)
```

    ##       estimator    param        Est         SE   bootstrap bootstraps
    ## beta         ML     beta 0.13279939 0.13728035 multinomial         19
    ## gamma        ML     tau2 0.08112823 0.08448746 multinomial         19
    ## zeta1        ML lambda_1 0.54845336 0.61595727 multinomial         19
    ##       percentile_lower percentile_upper
    ## beta      -0.032735627        0.3765500
    ## gamma      0.001500518        0.1984563
    ## zeta1      0.088833781        3.0549676

## Related Work

We want to recognize other packages that provide functions to selection
modeling.

Few packages are available to estimate selection models assuming that
effects are independent. The `metafor` package now includes the
`selmodel()` function which allows users to fit different types of
selection models (Viechtbauer 2010). The `weightr` package includes
functions to estimate weight-function models described in Vevea and
Hedges (1995; Coburn and Vevea 2019). However, the functions available
in these packages can only be applied to meta-analytic data assuming
that the effects are independent.

The `PublicationBias` package provides sensitivity analyses for
publication bias that incorporates robust variance estimation
(Braginsky, Mathur, and VanderWeele 2023). The analyses allow for
dependent effect sizes in meta-analytic data. However, the approach that
is implemented in the package is based on pre-specified degree of
selective reporting and only works for the simplest form of the
step-function model.

## Acknowledgements

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-PublicationBias" class="csl-entry">

Braginsky, Mika, Maya Mathur, and Tyler J. VanderWeele. 2023.
*PublicationBias: Sensitivity Analysis for Publication Bias in
Meta-Analyses*. <https://CRAN.R-project.org/package=PublicationBias>.

</div>

<div id="ref-weightr" class="csl-entry">

Coburn, Kathleen M., and Jack L. Vevea. 2019. *Weightr: Estimating
Weight-Function Models for Publication Bias*.
<https://CRAN.R-project.org/package=weightr>.

</div>

<div id="ref-vevea1995" class="csl-entry">

Vevea, Jack L., and Larry V. Hedges. 1995. “A General Linear Model for
Estimating Effect Size in the Presence of Publication Bias.”
*Psychometrika* 60 (3): 419435. <https://doi.org/10.1007/BF02294384>.

</div>

<div id="ref-metafor" class="csl-entry">

Viechtbauer, Wolfgang. 2010. “Conducting Meta-Analyses in R with the
<span class="nocase">metafor</span> Package.” *Journal of Statistical
Software* 36 (3): 1–48. <https://doi.org/10.18637/jss.v036.i03>.

</div>

<div id="ref-metadat" class="csl-entry">

White, Thomas, Daniel Noble, Alistair Senior, W. Kyle Hamilton, and
Wolfgang Viechtbauer. 2022. *Metadat: Meta-Analysis Datasets*.
<https://CRAN.R-project.org/package=metadat>.

</div>

</div>
