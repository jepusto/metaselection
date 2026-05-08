# Plot the selection weights implied by an estimated selection model.

For a fitted model of class `"selmodel"`, create a plot of the selection
weights implied by the model parameter estimates. If the model includes
bootstrapped confidence intervals, then the plot will also display the
selection weights implied by each bootstrap replicate of the parameter
estimates.

## Usage

``` r
selection_plot(
  mod,
  limits = c(0, 1),
  pts = 200L,
  ref_pval = NULL,
  transform = "identity",
  expand = ggplot2::expansion(0, 0.01),
  ...
)

# S3 method for class 'selmodel'
selection_plot(
  mod,
  limits = c(0, 1),
  pts = 200L,
  ref_pval = NULL,
  transform = "identity",
  expand = ggplot2::expansion(0, 0.01),
  fill = "blue",
  alpha = 0.5,
  step_linetype = "dashed",
  ...
)

# S3 method for class 'boot.selmodel'
selection_plot(
  mod,
  limits = c(0, 1),
  pts = 200L,
  ref_pval = NULL,
  transform = "identity",
  expand = ggplot2::expansion(0, 0.01),
  color = "black",
  linewidth = 1.2,
  step_linetype = "dashed",
  draw_boots = TRUE,
  fill = "blue",
  alpha = 0.5,
  boot_color = "blue",
  boot_alpha = 0.1,
  ...
)
```

## Arguments

- mod:

  Fitted model of class `"selmodel"`.

- limits:

  numeric vector of length 2 specifying the minimum and maximum p-values
  to plot.

- pts:

  Number of points for which to calculate selection weights, with a
  default of 200 points, evenly spaced between the specified limits.

- ref_pval:

  Numeric value of a p-value at which to standardize the weights. If not
  `NULL`, then a p-value of `ref_pval` will have selection weight of 1
  and selection weights for all other p-values will be calculated
  relative to `ref_pval`.

- transform:

  Character string specifying the name of a transformation function or
  the transformation function itself, as defined in the scales package.
  The transform is passed to
  [`ggplot2::scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).
  The default transform is `"identity"`. Other useful transforms for
  p-values are `"sqrt"` for square-root or `"asn"` for the arc-sin
  square root.

- expand:

  Passed to the `expand` argument of
  [`ggplot2::scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).

- ...:

  further arguments passed to
  [`ggplot2::scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).

- fill:

  character string specifying the fill-color to use when `mod` does not
  include bootstrap replications, with a default of `"blue"`. Passed to
  [`ggplot2::geom_area()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).

- alpha:

  numeric value specifying the opacity of the filled area plot, with a
  default of 0.5. Passed to
  [`ggplot2::geom_area()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).
  Only used when `mod` does not include bootstrap replications.

- step_linetype:

  character string specifying the type of line to draw to indicate
  p-value thresholds assumed in `mod`.

- color:

  character string specifying the line color to use for drawing the
  estimated selection weights, with a default of `"black"`. Passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  Only used when `mod` includes bootstrap replications.

- linewidth:

  numeric value specifying the line width to use for drawing the
  estimated selection weights, with a default of 1.2. Passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  Only used when `mod` includes bootstrap replications.

- draw_boots:

  logical value indicating whether to draw the selection weights for
  each bootstrap replication, with a default of `TRUE`.

- boot_color:

  character string specifying the line color to use for drawing the
  selection weights of each bootstrap replication, with a default of
  `"blue"`. Passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  Only used when `mod` includes bootstrap replications.

- boot_alpha:

  numeric value specifying the opacity of the lines for drawing the
  selection weights of each bootstrap replication, with a default of
  `"blue"`. Passed to
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html).
  Only used when `mod` includes bootstrap replications.

## Value

A `ggplot2` object.

## Examples

``` r
mod <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = c(0.025, .5),
  estimator = "CML",
  bootstrap = "none"
)

selection_plot(mod, fill = "purple")


# rescale the horizontal axis using arc-sin square root
selection_plot(mod, fill = "purple", transform = "asn") 



mod_boot <- selection_model(
  data = self_control,
  yi = g,
  sei = se_g,
  cluster = studyid,
  steps = c(0.025, .5),
  estimator = "ARGL",
  bootstrap = "multinomial",
  CI_type = "percentile",
  R = 9
)

 selection_plot(mod_boot, transform = "sqrt")

 selection_plot(mod_boot, transform = "sqrt", draw_boots = FALSE) # turn off bootstrap lines

 selection_plot(mod_boot, transform = "sqrt", color = "red", boot_color = "orange") # change colors

```
