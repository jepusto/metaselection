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
#'   estimator = "ML",
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
#'   estimator = "ML",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 9
#' )
#'
#' selection_wts(mod_boot, pvals = seq(0, 1, 0.2))
#' 
#' 


selection_wts <- function(mod, pvals, ...) UseMethod("selection_wts")


#' @export

selection_wts.default <- function(mod, pvals, params, steps, ...) {
  
  mod <- match.arg(mod, c("step","beta"))
  
  if (min(pvals) < 0 | max(pvals) > 1) stop("pvals must be in the interval [0,1].")
  
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
  
  sel_fun(pvals)
}

#' @rdname selection_wts
#' @param bootstraps If \code{mod} includes bootstrap replications, then setting
#'   \code{bootstraps = TRUE} will return selection weights for each bootstrap
#'   replication, in addition to the selection weights implied by the model
#'   parameter estimates. Ignored if \code{mod} does not include bootstrap
#'   replications.
#'
#' @export

selection_wts.step.selmodel <- function(mod, pvals, bootstraps = TRUE, ...) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_wts() is not available for models that include moderators of the selection parameters.")
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)
  
  wts <- data.frame(
    p = pvals, 
    wt = sel_fun(pvals)
  )
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)

  R <- eval(mod$cl$R, envir = parent.frame())
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- stats::reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  boot_wts <- apply(sel_param_boots, 1, \(zeta) step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)(pvals), simplify = FALSE)
  
  boot_wts <- data.frame(
    rep = rep(1:R, each = length(pvals)),
    p = rep(pvals, times = R),
    wt = unlist(boot_wts)
  )
  
  return(list(wts = wts, boot_wts = boot_wts))
  
}

#' @rdname selection_wts
#' @export

selection_wts.beta.selmodel <- function(mod, pvals, bootstraps = TRUE, ...) {
  
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
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)
  
  R <- eval(mod$cl$R, envir = parent.frame())
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- stats::reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  boot_wts <- apply(
    sel_param_boots, 1, 
    \(zeta) beta_fun(
      delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
      trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
      renormalize = FALSE
    )(pvals), 
    simplify = FALSE
  )
  
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
#' @param pts Number of points for which to calculate selection weights, with a
#'   default of 200 points, evenly spaced between 0 and 1.
#' @param ... further arguments passed to some methods.
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
#'   estimator = "ML",
#'   bootstrap = "none"
#' )
#' 
#' selection_plot(mod, fill = "purple")
#' 
#' 
#' mod_boot <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   steps = c(0.025, .5),
#'   estimator = "ML",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 9
#' )
#' 
#'  selection_plot(mod_boot)
#'  selection_plot(mod_boot, draw_boots = FALSE) # turn off bootstrap lines
#'  selection_plot(mod_boot, color = "red", boot_color = "orange") # change colors
#' 

selection_plot <- function(mod, pts = 200L, ...) UseMethod("selection_plot")


#' @export

selection_plot.default <- function(mod, pts = 200L, ...) {
  mod_class <- paste(class(mod), collapse = ", ")
  msg <- paste0("There is no `selection_plot` method available for objects of class ", mod_class, ".")
  stop(msg)
}

#' @rdname selection_plot
#' @param fill character string specifying the fill-color to use when \code{mod}
#'   does not include bootstrap replications, with a default of \code{"blue"}.
#'   Passed to \code{ggplot2::geom_area()}.
#' @param alpha numeric value specifying the opacity of the filled area plot,
#'   with a default of 0.5. Passed to \code{ggplot2::geom_area()}. Only used
#'   when \code{mod} does not include bootstrap replications.
#' @param step_linetype character string specifying the type of line to draw to
#'   indicate p-value thresholds assumed in \code{mod}.
#' @export
#'
#' @importFrom rlang .data

selection_plot.selmodel <- function(
  mod, 
  pts = 200L, 
  fill = "blue",
  alpha = 0.5,
  step_linetype = "dashed",
  ...
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)
  steps <- mod$steps
  dat <- selection_wts(mod, pvals = pts, bootstrap = FALSE)
  
  ggplot2::ggplot(dat) + 
    ggplot2::aes(x = .data$p, y = .data$wt) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0, 0.01)) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,0))) + 
    ggplot2::geom_vline(xintercept = steps, linetype = step_linetype) + 
    ggplot2::geom_area(fill = fill, alpha = alpha) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
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
    pts = 200L, 
    color = "black",
    linewidth = 1.2, 
    step_linetype = "dashed",
    draw_boots = TRUE,
    boot_color = "blue",
    boot_alpha = 0.1,
    ...
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)
  steps <- mod$steps
  
  dat <- selection_wts(mod, pvals = pts, bootstraps = TRUE)
  R <- eval(mod$cl$R, envir = parent.frame())
  
  p <- ggplot2::ggplot(dat$wts) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0, 0.01)) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,NA))) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
    ) + 
    ggplot2::geom_vline(xintercept = steps, linetype = step_linetype) + 
    ggplot2::theme_minimal()
  
  if (draw_boots) {
    p <- p + ggplot2::geom_line(data = dat$boot_wts, ggplot2::aes(x = .data$p, y = .data$wt, group = .data$rep), color = boot_color, alpha = 100 * boot_alpha / R)
  }
  
  p + ggplot2::geom_line(ggplot2::aes(x = .data$p, y = .data$wt), color = color, linewidth = linewidth)
    
}
