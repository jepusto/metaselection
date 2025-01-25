
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
<!-- This can take place at the study level, with entire studies left unpublished, or at the effect size level, with only positive effects reported and those with negative or non-affirming effects omitted from the publication. -->
Selective reporting is a major concern for research syntheses because it
distorts the evidence base available for a meta-analysis, skewing
meta-analytic averages toward favorable findings and misrepresenting the
true population of effects. Failure to account for selective reporting
can lead to inflated effect size estimates from meta-analysis and biased
estimates of heterogeneity, making it difficult to draw accurate
conclusions from a synthesis.

There are many tools available already to investigate and correct for
selective reporting. Widely used methods include: graphical diagnostics
like funnel plots, tests and adjustments for funnel plot asymmetry like
trim-and-fill, Egger’s regression, PET/PEESE, selection models, and
$p$-value diagnostics. However, very few methods for investigating
selective reporting can accommodate dependent effect sizes. This
limitation poses a problem for meta-analyses in education, psychology
and other social sciences, where dependent effects are a common feature
of meta-analytic data.

Dependent effect sizes occur when primary studies report results for
multiple measures of an outcome construct, collect repeated measures of
an outcome across multiple time-points, or involve comparisons between
multiple intervention conditions. Ignoring the dependency of effect size
estimates included in a meta-analysis leads to overly narrow confidence
intervals, hypothesis tests with inflated type one error rates, and
incorrect inferences.
<!-- Numerous methods have been developed to account for effect size dependencies, and some of these have been combined with a few of the available techniques for investigating selective reporting. -->
<!-- However, these combined methods are currently limited to techniques based on regression adjustment or sensitivity analyses based on simple forms of selection. -->
X (2024) developed and examined methods for investigating and accounting
for selective reporting in meta-analytic models that also account for
dependent effect sizes. Their simulation results show that combining
selection models with robust variance estimation to account for
dependent effects reduces bias in the estimate of the overall effect
size. Combining the selection models with cluster bootstrapping leads to
confidence intervals with close-to-nominal coverage rates.

The metaselection package provides an implementation of several
meta-analytic selection models. The main function, `selection_model()`,
fits step function and beta density selection models. To handle
dependence in the effect size estimates, the function provides options
to use cluster-robust (sandwich) variance estimation or cluster
bootstrapping to assess uncertainty in the model parameter estimates.

## Installation

You can install the development version of the package, along with a
vignette demonstrating how to use it, from GitHub with:

``` r
# If not already installed, first run install.packages("remotes")

remotes::install_github("jepusto/metaselection", build_vignettes = TRUE)
```

It may take a few minutes to install the package and vignette. Setting
`build_vignettes = FALSE` will lead to faster installation, although it
will preclude viewing the package vignette.

## Example

The following example uses data from a meta-analysis by Lehmann
meta-analysis which examined the effects of color red on attractiveness
judgments. The dataset is included in the `metadat` package (White et
al. 2022) as `dat.lehmann`. In the code below, we fit a step function
selection model to the Lehmann dataset using the `selection_model()`
function, with confidence intervals computed using cluster
bootstrapping. For further details, please see the vignette.

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
  # Set R to a much higher number of bootstrap replications, 
  # such as 1999, to obtain confidence intervals with 
  # more accurate coverage rates
)

summary(mod_3PSM_boot)
```

    ## Step Function Model with Cluster Bootstrapping 
    ##  
    ## Call: 
    ## selection_model(data = dat.lehmann2018, yi = yi, sei = sei, cluster = study, 
    ##     selection_type = "step", steps = 0.025, CI_type = "percentile", 
    ##     bootstrap = "multinomial", R = 19)
    ## 
    ## Number of clusters = 41; Number of effects = 81
    ## 
    ## Steps: 0.025 
    ## Estimator: composite marginal likelihood 
    ## Variance estimator: robust 
    ## Bootstrap type: multinomial 
    ## Number of bootstrap replications: 19 
    ## 
    ## Log composite likelihood of selection model: -44.46436
    ## Inverse selection weighted partial log likelihood: 58.35719 
    ## 
    ## Mean effect estimates:                                               
    ##                            Percentile Bootstrap
    ##  Coef. Estimate Std. Error      Lower     Upper
    ##   beta    0.133      0.137    -0.0327     0.377
    ## 
    ## Heterogeneity estimates:                                               
    ##                            Percentile Bootstrap
    ##  Coef. Estimate Std. Error      Lower     Upper
    ##   tau2   0.0811     0.0845     0.0015     0.198
    ## 
    ## Selection process estimates:
    ##  Step: 0 < p <= 0.025; Studies: 16; Effects: 25                                                 
    ##                              Percentile Bootstrap
    ##    Coef. Estimate Std. Error      Lower     Upper
    ##  lambda0        1        ---        ---       ---
    ## 
    ##  Step: 0.025 < p <= 1; Studies: 29; Effects: 56                                                 
    ##                              Percentile Bootstrap
    ##    Coef. Estimate Std. Error      Lower     Upper
    ##  lambda1    0.548      0.616     0.0888      3.05

The beta estimate of 0.133, with a 95% confidence interval -0.033,
0.377, represents the overall average effect after accounting for both
selection bias and dependent effects. Tau (0.081) is the estimated total
variance, including both between- and within-study heterogeneity.
`lambda1` is the selection parameter. The estimate of 0.548 indicates
that effect size estimates with one-sided $p$-values greater than 0.025
are only about half as likely to be reported as estimates that are
positive and statistically significant (i.e., estimates with
$p < 0.025$).

The package is designed to work with the `progressr` package. To turn on
progress bars for all bootstrap calculations, use

``` r
progressr::handlers(global = TRUE)
```

See `vignette("progressr-intro")` for further details.

The package is also designed to work with the `future` package for
parallel computing. To enable parallel computation of bootstrap
calculations, simply set an appropriate parallelization plan such as

``` r
library(future)
plan(multisession)
```

The `metaselection` package vignette includes a more detailed
demonstration.

## Related Work

Several existing meta-analysis packages provide implementations of
selection models, but are limited by the assumption of independent
effects. The `metafor` package (Viechtbauer 2010) includes the
`selmodel()` function, which allows users to fit many different types of
selection models. The `weightr` package (Coburn and Vevea 2019) includes
functions to estimate a class of $p$-value selection models described in
Vevea and Hedges (1995). However, because these packages assume effect
sizes to be independent, the results they produce will have incorrect
standard errors and misleadingly narrow confidence intervals for
datasets containing multiple effects drawn from the same study. In
addition, the `PublicationBias` package (Braginsky, Mathur, and
VanderWeele 2023) implements sensitivity analyses for selective
reporting bias that incorporate cluster-robust variance estimation
methods for handling dependent effect sizes. However, the sensitivity
analyses implemented in the package are based on a pre-specified degree
of selective reporting, rather than allowing the degree of selection to
be estimated from the data. The sensitivity analyses are also based on a
specific and simple form of selection model, and do not allow
consideration of more complex forms of selection functions. The
`metaselection` package goes beyond these other tools both by
considering more complex forms of selective reporting and by correcting
for selective reporting bias while accommodating meta-analytic datasets
that include dependent effect sizes.

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

<div id="ref-vevea1995general" class="csl-entry">

Vevea, Jack L, and Larry V Hedges. 1995. “A General Linear Model for
Estimating Effect Size in the Presence of Publication Bias.”
*Psychometrika* 60 (3): 419–35. <https://doi.org/10.1007/BF02294384>.

</div>

<div id="ref-Viechtbauer2010conducting" class="csl-entry">

Viechtbauer, Wolfgang. 2010. “<span class="nocase">Conducting
meta-analyses in R with the metafor package</span>.” *Journal of
Statistical Software* 36 (3): 1–48.
<https://doi.org/10.18637/jss.v036.i03>.

</div>

<div id="ref-metadat" class="csl-entry">

White, Thomas, Daniel Noble, Alistair Senior, W. Kyle Hamilton, and
Wolfgang Viechtbauer. 2022. *Metadat: Meta-Analysis Datasets*.
<https://CRAN.R-project.org/package=metadat>.

</div>

</div>
