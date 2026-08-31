
<!-- badges: start -->

[![R-CMD-check](https://github.com/jepusto/metaselection/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jepusto/metaselection/actions/workflows/R-CMD-check.yaml)
[![Codecov
Status](https://codecov.io/gh/jepusto/metaselection/graph/badge.svg?token=8T7IUFT1QV)](https://codecov.io/gh/jepusto/metaselection)
[![CRAN
Version](http://www.r-pkg.org/badges/version/metaselection)](https://CRAN.R-project.org/package=metaselection)
[![](http://cranlogs.r-pkg.org/badges/grand-total/metaselection)](https://CRAN.R-project.org/package=metaselection)
[![](http://cranlogs.r-pkg.org/badges/last-month/metaselection)](https://CRAN.R-project.org/package=metaselection)

<!-- badges: end -->

# metaselection

Selective reporting is a major concern for research syntheses because it
distorts the evidence base available for a meta-analysis. There are many
tools available to investigate and correct for selective reporting but
few of them can accommodate dependent effect sizes. Ignoring the
dependency of effect size estimates included in a meta-analysis leads to
overly narrow confidence intervals, hypothesis tests with inflated Type
1 error rates, and incorrect inferences.

[Pustejovsky et al. (2025)](https://osf.io/preprints/metaarxiv/qg5x6_v4)
and [Citkowicz et al.
(2026)](https://osf.io/preprints/metaarxiv/wjpxk_v1) developed methods
for investigating and accounting for selective reporting in
meta-analytic models that also account for dependent effect sizes. Their
simulation results show that using a marginal selection model reduces
bias in the estimate of the overall effect size compared to using
conventional summary meta-analysis models or techniques such as the
PET/PEESE adjustment. Moreover, combining the marginal selection model
with cluster bootstrapping — particularly the two-stage cluster
bootstrapping — leads to confidence intervals with close-to-nominal
coverage rates.

The `metaselection` package provides an implementation of the methods
proposed and evaluated by Pustejovsky et al. (2025) and Citkowicz et al.
(2026). The main function,
[`selection_model()`](https://jepusto.github.io/metaselection/reference/selection_model.html),
can fit step- and beta-function selection models. To handle dependence
in the effect size estimates, the function provides options to use
cluster-robust (sandwich) variance estimation or cluster bootstrapping
to assess uncertainty in the model parameter estimates. We highly
recommend using cluster bootstrapping, particularly the two-stage
cluster bootstrapping with percentile bootstrap confidence intervals.

## Installation

You can install the latest release of the package from CRAN by running:

``` r
install.packages("metaselection")
```

You can install the development version of the package, along with a
vignette demonstrating how to use it, from GitHub with:

``` r
# If not already installed, first run install.packages("remotes")

remotes::install_github("jepusto/metaselection", dependencies = TRUE, build_vignettes = TRUE)
```

It may take a few minutes to install the package and vignette. Setting
`build_vignettes = FALSE` will lead to faster installation, although it
will preclude viewing the package vignette.

After installation, you can view the vignette by running the command:

``` r
vignette("selection-models", package = "metaselection")
```

You can also [view the vignette
online](https://jepusto.github.io/metaselection/articles/selection-models.html)
on the `metaselection` package website.

## Example

The following example uses data from a meta-analysis by Lehmann et al.
(2018) which examined the effects of color red on attractiveness
judgments. The dataset is included in the `metadat` package (White et
al. 2022) as `dat.lehmann`. In the code below, we fit a three-parameter
step-function selection model to the Lehmann dataset using the
`selection_model()` function, with confidence intervals computed using
two-stage cluster bootstrapping. Additionally, we use a default set of
weak priors to regularize the estimates (for detail, see
[`define_priors()`](https://jepusto.github.io/metaselection/reference/define_priors.html)).
The function also has an argument that can be used to set the direction
of the alternative hypothesis used in computing the p-values (default
`alternative = "greater"`). For further details on how to use the
`selection_model()` function, please see the
[vignette](https://jepusto.github.io/metaselection/articles/selection-models.html).

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
  priors = define_priors(),
  alternative = "greater",
  CI_type = "percentile",
  bootstrap = "two-stage",
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
    ##     selection_type = "step", alternative = "greater", steps = 0.025, 
    ##     priors = define_priors(), CI_type = "percentile", bootstrap = "two-stage", 
    ##     R = 19)
    ## 
    ## Number of clusters = 41; Number of effects = 81
    ## 
    ## Steps: 0.025 
    ## Estimator: composite marginal likelihood 
    ## Variance estimator: robust 
    ## Bootstrap type: two-stage 
    ## Number of bootstrap replications: 19 
    ## 
    ## Log composite likelihood of selection model: -44.46655
    ## Inverse selection weighted partial log likelihood: 59.53697 
    ## 
    ## Mean effect estimates:                                               
    ##                            Percentile Bootstrap
    ##  Coef. Estimate Std. Error      Lower     Upper
    ##   beta    0.131      0.135    -0.0492     0.412
    ## 
    ## Heterogeneity estimates:                                               
    ##                            Percentile Bootstrap
    ##  Coef. Estimate Std. Error      Lower     Upper
    ##   tau2   0.0794     0.0815    0.00298     0.223
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
    ##  lambda1     0.54      0.601     0.0844      4.35

The beta estimate of 0.131, with a 95% confidence interval of \[-0.049,
0.412\], represents the overall average effect after accounting for both
selection bias and dependent effects. The `tau2` ($\tau^2$) estimate of
0.079 is the estimated total variance, including both between- and
within-study heterogeneity. `lambda1` ($\lambda_1$) is the selection
parameter. The estimate of 0.54 indicates that effect size estimates
with one-sided $p$-values greater than 0.025 are only about half as
likely to be reported as estimates that are positive and statistically
significant (i.e., estimates with $p < 0.025$). [This
infographic](https://www.air.org/sites/default/files/2025-09/How-to-Read-Step-Function-Selection-Model-Results-infographic-Sept-2025.pdf)
provides further guidance on interpreting the model output.
<!--# QUESTION: Do you want to move this infographic to our package website as a separate vignette so we don't have to rely on AIR for keeping the page active? QUESTION: do you know why the results are slightly different from the numbers in the current website readme?  -->

## Parallel computing and tracking progress

Bootstrapping can be computationally intensive. To reduce computation
time, we designed this package to work with the
[`future`](https://future.futureverse.org/) package for parallel
computing (Bengtsson 2021). To enable parallel computation of bootstrap
calculations, simply set an appropriate parallelization plan such as

``` r
library(future)
plan(multisession)
```

Our
[vignette](https://jepusto.github.io/metaselection/articles/selection-models.html#parallel-processing)
includes a more detailed demonstration.

The package is also designed to work with the
[`progressr`](https://progressr.futureverse.org/) package (Bengtsson
2026). To turn on progress bars for all bootstrap calculations, use

``` r
progressr::handlers(global = TRUE)
```

See
[`vignette("progressr-intro")`](https://progressr.futureverse.org/articles/progressr-01-intro.html)
for further details.

## Related Work

Several existing meta-analysis packages provide implementations of
selection models, but are limited by the assumption of independent
effects:

- The `metafor` package (Viechtbauer 2010) includes the `selmodel()`
  function, which allows users to fit many different types of selection
  models.
- The `weightr` package (Coburn and Vevea 2019) includes functions to
  estimate a class of $p$-value selection models described in Vevea and
  Hedges (1995). However, because these packages assume effect sizes to
  be independent, the results they produce will have incorrect standard
  errors and misleadingly narrow confidence intervals for datasets
  containing multiple effects drawn from the same study.
- The `RoBMA` package (Bartoš and Maier 2020) implements Bayesian
  ensemble models that include step-function selection models as one
  component of the ensemble. `RoBMA` includes an implementation of a
  multilevel extension to the step-function selection model, as
  described in Bartoš et al. (2026), a different approach for handling
  dependence that is an alternative to the marginal selection models
  implemented in `metaselection`.
- The `PublicationBias` package (Braginsky et al. 2023) implements
  sensitivity analyses for selective reporting bias that incorporate
  cluster-robust variance estimation methods for handling dependent
  effect sizes. However, the sensitivity analyses implemented in the
  package are based on a pre-specified degree of selective reporting,
  rather than allowing the degree of selection to be estimated from the
  data. The sensitivity analyses are also based on a specific and simple
  form of selection model, and do not allow consideration of more
  complex forms of selection functions.

The `metaselection` package goes beyond these other tools both by
considering more complex forms of selective reporting and by correcting
for selective reporting bias while accommodating meta-analytic datasets
that include dependent effect sizes.

## Acknowledgements

The development of this software was supported, in whole or in part, by
the Institute of Education Sciences, U.S. Department of Education,
through grant
[R305D220026](https://ies.ed.gov/funding/grantsearch/details.asp?ID=5730)
to the American Institutes for Research. The opinions expressed are
those of the authors and do not represent the views of the Institute or
the U.S. Department of Education.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-RoBMA" class="csl-entry">

Bartoš, František, and Maximilian Maier. 2020. *RoBMA: An r Package for
Robust Bayesian Meta-Analyses*.
<https://CRAN.R-project.org/package=RoBMA>.

</div>

<div id="ref-bartovs2026robust" class="csl-entry">

Bartoš, František, Maximilian Maier, and Eric-Jan Wagenmakers. 2026.
“Robust Bayesian Multilevel Meta-Analysis: Adjusting for Publication
Bias in the Presence of Dependent Effect Sizes.” *Behavior Research
Methods* 58 (6): 165.

</div>

<div id="ref-future" class="csl-entry">

Bengtsson, Henrik. 2021. “A Unifying Framework for Parallel and
Distributed Processing in r Using Futures.” *The R Journal* 13 (2):
208–27. <https://doi.org/10.32614/RJ-2021-048>.

</div>

<div id="ref-progressr" class="csl-entry">

Bengtsson, Henrik. 2026. *Progressr: An Inclusive, Unifying API for
Progress Updates*. <https://doi.org/10.32614/CRAN.package.progressr>.

</div>

<div id="ref-PublicationBias" class="csl-entry">

Braginsky, Mika, Maya Mathur, and Tyler J. VanderWeele. 2023.
*PublicationBias: Sensitivity Analysis for Publication Bias in
Meta-Analyses*. <https://CRAN.R-project.org/package=PublicationBias>.

</div>

<div id="ref-citkowicz2026estimating" class="csl-entry">

Citkowicz, Martyna, James E. Pustejovsky, and Megha Joshi. 2026.
*Estimating Beta-Function Selection Models in Meta-Analysis with
Dependent Effects*. <https://doi.org/10.31222/osf.io/wjpxk_v1>.

</div>

<div id="ref-weightr" class="csl-entry">

Coburn, Kathleen M., and Jack L. Vevea. 2019. *Weightr: Estimating
Weight-Function Models for Publication Bias*.
<https://CRAN.R-project.org/package=weightr>.

</div>

<div id="ref-lehmann2018meta" class="csl-entry">

Lehmann, Gabrielle K, Andrew J Elliot, and Robert J Calin-Jageman. 2018.
“Meta-Analysis of the Effect of Red on Perceived Attractiveness.”
*Evolutionary Psychology* 16 (4): 1474704918802412.

</div>

<div id="ref-pustejovsky2025estimation" class="csl-entry">

Pustejovsky, James E., Martyna Citkowicz, and Megha Joshi. 2025.
*Estimation and Inference for Step-Function Selection Models in
Meta-Analysis with Dependent Effects*.
<https://doi.org/10.31222/osf.io/qg5x6_v1>.

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
