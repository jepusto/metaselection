# Summarize results from a `selmodel` object

Summarize relevant results from a `selmodel` object.

## Usage

``` r
# S3 method for class 'selmodel'
summary(object, transf_gamma = TRUE, transf_zeta = TRUE, digits = 3, ...)
```

## Arguments

- object:

  Fitted model of class `"selmodel"`.

- transf_gamma:

  logical with `TRUE` (the default) indicating that the heterogeneity
  parameter estimates (called gamma) should be transformed by
  exponentiating.

- transf_zeta:

  logical with `TRUE` (the default) indicating that the selection
  parameter estimates (called zeta) should be transformed by
  exponentiating.

- digits:

  Minimum number of significant digits to be used, with a default of 3.

- ...:

  further arguments passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

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
summary(res_ML, transf_gamma = FALSE, transf_zeta = FALSE)
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
#>                                    Large Sample
#>  Coef. Estimate Std. Error p-value Lower  Upper
#>  gamma    -3.23      0.725     --- -4.66  -1.81
#> 
#> Selection process estimates:
#>  Step: 0 < p <= 0.025; Studies: 18; Effects: 32                                               
#>                                    Large Sample
#>  Coef. Estimate Std. Error p-value Lower  Upper
#>  zeta0        0        ---     ---   ---    ---
#> 
#>  Step: 0.025 < p <= 1; Studies: 24; Effects: 126                                                
#>                                     Large Sample
#>  Coef. Estimate Std. Error p-value  Lower  Upper
#>  zeta1   0.0344      0.489   0.944 -0.924  0.992
```
