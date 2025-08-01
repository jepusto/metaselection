#-------------------------------------------------------------------------------
# Calculate selection weights for specified p-values

#' @title Calculate model-implied weights for specified p-values.
#'
#' @description Calculates the selection weights implied by an estimated model
#'   for a user-specified p-value or set of p-values.
#'
#' @param mod Fitted model of class \code{"selmodel"}.
#' @param pvals Numeric vector of p-values for which to calculate selection
#'   weights.
#' @param ref_pval Numeric value of a p-value at which to standardize the
#'   weights. If not \code{NULL}, then a p-value of \code{ref_pval} will have
#'   selection weight of 1 and selection weights for all other p-values will be
#'   calculated relative to \code{ref_pval}.
#' @param ... further arguments passed to some methods.
#'
#' @returns If \code{mod} does not include bootstrapped confidence intervals or
#'   if the argument \code{bootstraps = FALSE}, then \code{selection_wts} will
#'   return a \code{data.frame} containing the user-specified p-values and the
#'   selection weights implied by the estimated model parameters.
#'
#'   If \code{mod} does include bootstrapped confidence intervals (i.e., when
#'   \code{inherits(mod, "boot.selmodel")} is \code{TRUE}) and the argument
#'   \code{bootstraps = TRUE}, then \code{selection_wts} will return a list with
#'   two elements. The first element is a \code{data.frame} containing the
#'   user-specified p-values and the selection weights implied by the estimated
#'   model parameters. The second element is a \code{data.frame} containing the
#'   user-specified p-values and the selection weights implied by each bootstrap
#'   replicate of the model parameter estimates. The \code{data.frame} includes
#'   an additional variable, \code{rep}, identifying the bootstrap replicate.
#'
#' @export
#'
#' @examples
#' mod <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = c(0.025, .5),
#'   estimator = "CML",
#'   bootstrap = "none"
#' )
#'
#' selection_wts(mod, pvals = seq(0, 1, 0.2))
#'
#' mod_boot <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = c(0.025, .5),
#'   estimator = "CML",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 9
#' )
#'
#' selection_wts(mod_boot, pvals = seq(0, 1, 0.2))
#'
#' 


selection_wts <- function(mod, pvals, ref_pval, ...) UseMethod("selection_wts")


#' @export

selection_wts.default <- function(mod, pvals, ref_pval, params, steps, ...) {
  
  mod <- match.arg(mod, c("step","beta"))
  
  if (min(pvals) < 0 | max(pvals) > 1) stop("`pvals` must be in the interval [0,1].")
  if (length(ref_pval) > 1) stop("`ref_pval` must be a single value.")
  if (ref_pval < 0 | ref_pval > 1) stop("`ref_pval` must be between 0 and 1.")
  
  if (mod == "step") {
    
    if (length(params) != length(steps)) stop("The length of `params` must be equal to the length of `steps`.")
    sel_fun <- step_fun(
      cut_vals = steps, 
      weights = exp(params), 
      renormalize = FALSE
    )
    
  } else if (mod == "beta") {
    
    if (length(params) != 2L) stop("The `params` argument must be a vector of length 2.")
    if (length(steps) != 2L) stop("The `steps` argument must be a vector of length 2.")
    
    sel_fun <- beta_fun(
      delta_1 = exp(params[1]), delta_2 = exp(params[2]),
      trunc_1 = steps[1], trunc_2 = steps[2],
      renormalize = FALSE
    )
  }
  
  if (missing(ref_pval)) {
    sel_fun(pvals)
  } else {
    sel_fun(pvals) / sel_fun(ref_pval)
  }
  
}

#' @rdname selection_wts
#' @param bootstraps If \code{mod} includes bootstrap replications, then setting
#'   \code{bootstraps = TRUE} will return selection weights for each bootstrap
#'   replication, in addition to the selection weights implied by the model
#'   parameter estimates. Ignored if \code{mod} does not include bootstrap
#'   replications.
#'
#' @export

selection_wts.step.selmodel <- function(mod, pvals = NULL, ref_pval = NULL, bootstraps = TRUE, ...) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_wts() is not available for models that include moderators of the selection parameters.")
  if (length(ref_pval) > 1) stop("`ref_pval` must be a single value.")
  
  if (is.null(pvals)) {
    if (is.null(mod$cl_pi)) {
      yi <- eval(mod$cl$yi, envir = mod$mf)
      sei <- eval(mod$cl$sei, envir = mod$mf)
      pvals <- pnorm(yi / sei, lower.tail = FALSE)
    } else {
      pvals <- eval(mod$cl$pi, envir = mod$mf)
    }
  }
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)
  
  wts <- data.frame(
    p = pvals, 
    wt = sel_fun(pvals)
  )
  
  if (!is.null(ref_pval)) {
    if (ref_pval < 0 | ref_pval > 1) stop("`ref_pval` must be between 0 and 1.")
    ref_wt <- sel_fun(ref_pval)
    wts$wt <- wts$wt / ref_wt
  }
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)

  R <- mod$R
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- stats::reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  if (is.null(ref_pval)) {
    boot_wts <- apply(
      sel_param_boots, 1, \(zeta) 
      step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)(pvals), 
      simplify = FALSE
    )
  } else {
    boot_wts <- apply(
      sel_param_boots, 1, \(zeta) 
      step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)(pvals) / step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)(ref_pval), 
      simplify = FALSE
    )
  }
  
  boot_wts <- data.frame(
    rep = rep(1:R, each = length(pvals)),
    p = rep(pvals, times = R),
    wt = unlist(boot_wts)
  )
  
  return(list(wts = wts, boot_wts = boot_wts))
  
}

#' @rdname selection_wts
#' @export

selection_wts.beta.selmodel <- function(mod, pvals = NULL, ref_pval = NULL, bootstraps = TRUE, ...) {
  
  if (length(ref_pval) > 1) stop("`ref_pval` must be a single value.")
  
  if (is.null(pvals)) {
    if (is.null(mod$cl_pi)) {
      yi <- eval(mod$cl$yi, envir = mod$mf)
      sei <- eval(mod$cl$sei, envir = mod$mf)
      pvals <- pnorm(yi / sei, lower.tail = FALSE)
    } else {
      pvals <- eval(mod$cl$pi, envir = mod$mf)
    }
  }
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- beta_fun(
    delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]),
    trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
    renormalize = FALSE
  )
  
  wts <- data.frame(
    p = pvals, 
    wt = sel_fun(pvals)
  )
  
  if (!is.null(ref_pval)) {
    if (ref_pval < 0 | ref_pval > 1) stop("`ref_pval` must be between 0 and 1.")
    ref_wt <- sel_fun(ref_pval)
    wts$wt <- wts$wt / ref_wt
  }
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)
  
  R <- mod$R
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- stats::reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  

  
  if (is.null(ref_pval)) {
    boot_wts <- apply(
      sel_param_boots, 1, 
      \(zeta) beta_fun(
        delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
        trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
        renormalize = FALSE
      )(pvals), 
      simplify = FALSE
    )
  } else {
    boot_wts <- apply(
      sel_param_boots, 1, 
      \(zeta) 
      beta_fun(
        delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
        trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
        renormalize = FALSE
      )(pvals) / 
      beta_fun(
        delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
        trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
        renormalize = FALSE
      )(ref_pval), 
      simplify = FALSE
    )
  }
  
  boot_wts <- data.frame(
    rep = rep(1:R, each = length(pvals)),
    p = rep(pvals, times = R),
    wt = unlist(boot_wts)
  )
  
  return(list(wts = wts, boot_wts = boot_wts))
}




#-------------------------------------------------------------------------------
# Plotting methods

#' @title Plot the selection weights implied by an estimated selection model.
#'
#' @description For a fitted model of class \code{"selmodel"}, create a plot of
#'   the selection weights implied by the model parameter estimates. If the
#'   model includes bootstrapped confidence intervals, then the plot will also
#'   display the selection weights implied by each bootstrap replicate of the
#'   parameter estimates.
#'
#' @param mod Fitted model of class \code{"selmodel"}.
#' @param limits numeric vector of length 2 specifying the minimum and maximum p-values to plot.
#' @param pts Number of points for which to calculate selection weights, with a
#'   default of 200 points, evenly spaced between the specified limits.
#' @param transform Character string specifying the name of a transformation function or the transformation function itself, as defined in the scales package. The transform is passed to \code{ggplot2::scale_x_continuous}. The default transform is \code{"identity"}. Other useful transforms for p-values are \code{"sqrt"} for square-root or \code{"asn"} for the arc-sin square root.
#' @param expand Passed to the \code{expand} argument of \code{ggplot2::scale_x_continuous}.
#' @param fill character string specifying the fill-color to use when \code{mod}
#'   does not include bootstrap replications, with a default of \code{"blue"}.
#'   Passed to \code{ggplot2::geom_area()}.
#' @param alpha numeric value specifying the opacity of the filled area plot,
#'   with a default of 0.5. Passed to \code{ggplot2::geom_area()}. Only used
#'   when \code{mod} does not include bootstrap replications.
#' @inheritParams selection_wts
#' @param ... further arguments passed to \code{ggplot2::scale_x_continuous}.
#'
#' @returns A \code{ggplot2} object.
#'
#' @export
#' 
#' @examples
#' mod <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = c(0.025, .5),
#'   estimator = "CML",
#'   bootstrap = "none"
#' )
#' 
#' selection_plot(mod, fill = "purple")
#' 
#' # rescale the horizontal axis using arc-sin square root
#' selection_plot(mod, fill = "purple", transform = "asn") 
#' 
#' 
#' mod_boot <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = c(0.025, .5),
#'   estimator = "ARGL",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 9
#' )
#' 
#'  selection_plot(mod_boot, transform = "sqrt")
#'  selection_plot(mod_boot, transform = "sqrt", draw_boots = FALSE) # turn off bootstrap lines
#'  selection_plot(mod_boot, transform = "sqrt", color = "red", boot_color = "orange") # change colors
#' 

selection_plot <- function(mod, limits = c(0,1), pts = 200L, ref_pval = NULL, transform = "identity", expand = ggplot2::expansion(0, 0.01), ...) UseMethod("selection_plot")


#' @export

selection_plot.default <- function(mod, limits = c(0,1), pts = 200L, ref_pval = NULL, transform = "identity", expand = ggplot2::expansion(0, 0.01), ...) {
  mod_class <- paste(class(mod), collapse = ", ")
  msg <- paste0("There is no `selection_plot` method available for objects of class ", mod_class, ".")
  stop(msg)
}

#' @rdname selection_plot
#' @param step_linetype character string specifying the type of line to draw to
#'   indicate p-value thresholds assumed in \code{mod}.
#' @export
#'
#' @importFrom rlang .data

selection_plot.selmodel <- function(
  mod, 
  limits = c(0,1),
  pts = 200L, 
  ref_pval = NULL,
  transform = "identity",
  expand = ggplot2::expansion(0, 0.01),
  fill = "blue",
  alpha = 0.5,
  step_linetype = "dashed",
  ...
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  if (identical(transform, "identity")) {
    pts <- seq(from = limits[1], to = limits[2], length.out = pts)
  } else {
    transf <- scales::as.transform(transform)
    limits_trans <- transf$transform(limits)
    pts_trans <- seq(from = limits_trans[1], to = limits_trans[2], length.out = pts)
    pts <- transf$inverse(pts_trans)
  }
  
  steps <- mod$steps
  dat <- selection_wts(mod, pvals = pts, ref_pval = ref_pval, bootstrap = FALSE)
  
  ggplot2::ggplot(dat) + 
    ggplot2::aes(x = .data$p, y = .data$wt) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = limits, expand = expand, transform = transform, ...) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,0))) + 
    ggplot2::geom_vline(xintercept = steps, linetype = step_linetype) + 
    ggplot2::geom_area(fill = fill, alpha = alpha) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Relative probability of selection"
    )
  
}

#' @rdname selection_plot
#' @param color character string specifying the line color to use for drawing
#'   the estimated selection weights, with a default of \code{"black"}. Passed
#'   to \code{ggplot2::geom_line()}. Only used when \code{mod} includes
#'   bootstrap replications.
#' @param linewidth numeric value specifying the line width to use for drawing
#'   the estimated selection weights, with a default of 1.2. Passed to
#'   \code{ggplot2::geom_line()}. Only used when \code{mod} includes bootstrap
#'   replications.
#' @param step_linetype character string specifying the type of line to draw to
#'   indicate p-value thresholds assumed in \code{mod}.
#' @param draw_boots logical value indicating whether to draw the selection
#'   weights for each bootstrap replication, with a default of \code{TRUE}.
#' @param boot_color character string specifying the line color to use for
#'   drawing the selection weights of each bootstrap replication, with a default
#'   of \code{"blue"}. Passed to \code{ggplot2::geom_line()}. Only used when
#'   \code{mod} includes bootstrap replications.
#' @param boot_alpha numeric value specifying the opacity of the lines for
#'   drawing the selection weights of each bootstrap replication, with a default
#'   of \code{"blue"}. Passed to \code{ggplot2::geom_line()}. Only used when
#'   \code{mod} includes bootstrap replications.
#'
#' @export

selection_plot.boot.selmodel <- function(
    mod, 
    limits = c(0,1),
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
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  if (identical(transform, "identity")) {
    pts <- seq(from = limits[1], to = limits[2], length.out = pts)
  } else {
    transf <- scales::as.transform(transform)
    limits_trans <- transf$transform(limits)
    pts_trans <- seq(from = limits_trans[1], to = limits_trans[2], length.out = pts)
    pts <- transf$inverse(pts_trans)
  }
  
  steps <- mod$steps
  
  dat <- selection_wts(mod, pvals = pts, ref_pval = ref_pval, bootstraps = TRUE)
  R <- mod$R
  
  p <- ggplot2::ggplot(dat$wts) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = limits, expand = expand, transform = transform, ...) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,NA))) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Relative probability of selection"
    ) + 
    ggplot2::geom_vline(xintercept = steps, linetype = step_linetype) + 
    ggplot2::theme_minimal()
  
  if (draw_boots) {
    p <- p + 
      ggplot2::geom_line(
        data = dat$boot_wts, 
        ggplot2::aes(x = .data$p, y = .data$wt, group = .data$rep), 
        color = boot_color, 
        alpha = 100 * boot_alpha / R
      ) +
      ggplot2::geom_line(
        ggplot2::aes(x = .data$p, y = .data$wt),
        color = color,
        linewidth = linewidth
      )
  } else {
    p <- p + 
      ggplot2::geom_area(
        ggplot2::aes(x = .data$p, y = .data$wt), 
        fill = fill, 
        alpha = alpha
      ) 
  }
  
  p
}
