#-------------------------------------------------------------------------------
# Calculate selection weights for specified p-values

#' @export

selection_wts <- function(mod, pts, ...) UseMethod("selection_wts")

#' @export

selection_wts.default <- function(mod, pts, params, steps) {
  
  mod <- match.arg(mod, c("step","beta"))
  
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
  
  sel_fun(pts)
}

#' @export

selection_wts.step.selmodel <- function(mod, pts, bootstraps = TRUE) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_wts() is not available for models that include moderators of the selection parameters.")
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)
  
  wts <- data.frame(
    p = pts, 
    wt = sel_fun(pts)
  )
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)

  R <- eval(mod$cl$R, envir = parent.frame())
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  boot_wts <- apply(sel_param_boots, 1, \(zeta) step_fun(cut_vals = mod$steps, weights = exp(zeta), renormalize = FALSE)(pts), simplify = FALSE)
  
  boot_wts <- data.frame(
    rep = rep(1:R, each = length(pts)),
    p = rep(pts, times = R),
    wt = unlist(boot_wts)
  )
  
  return(list(wts = wts, boot_wts = boot_wts))
  
}

#' @export

selection_wts.beta.selmodel <- function(mod, pts, bootstraps = TRUE) {
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- beta_fun(
    delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]),
    trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
    renormalize = FALSE
  )
  
  wts <- data.frame(
    p = pts, 
    wt = sel_fun(pts)
  )
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)
  
  R <- eval(mod$cl$R, envir = parent.frame())
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  boot_wts <- apply(
    sel_param_boots, 1, 
    \(zeta) beta_fun(
      delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
      trunc_1 = mod$steps[1], trunc_2 = mod$steps[2],
      renormalize = FALSE
    )(pts), 
    simplify = FALSE
  )
  
  boot_wts <- data.frame(
    rep = rep(1:R, each = length(pts)),
    p = rep(pts, times = R),
    wt = unlist(boot_wts)
  )
  
  return(list(wts = wts, boot_wts = boot_wts))
}




#-------------------------------------------------------------------------------
# Plotting methods

#' @export

selection_plot <- function(mod, pts, ...) UseMethod("selection_plot")

#' @export

selection_plot.default <- function(mod, pts, ...) {
  mod_class <- paste(class(mod), collapse = ", ")
  msg <- paste0("There is no `selection_plot` method available for objects of class ", mod_class, ".")
  stop(msg)
}

#' @export

selection_plot.selmodel <- function(
  mod, 
  pts = 200L, 
  color = "blue",
  alpha = 0.5,
  ...
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)

  dat <- selection_wts(mod, pts = pts, bootstrap = FALSE)
  
  ggplot2::ggplot(dat) + 
    ggplot2::aes(x = p, y = wt) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0, 0.01)) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,0))) + 
    ggplot2::geom_area(fill = color, alpha = alpha) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
    )
  
}

#' @export

selection_plot.boot.selmodel <- function(
    mod, 
    pts = 200L, 
    color = "black",
    linewidth = 1.2, 
    draw_boots = TRUE,
    boot_color = "blue",
    boot_alpha = 0.1,
    ...
) {
  
  if (!is.null(mod$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)
  
  dat <- selection_wts(mod, pts = pts, bootstraps = TRUE)
  R <- eval(mod$cl$R, envir = parent.frame())
  
  p <- ggplot2::ggplot(dat$wts) + 
    ggplot2::expand_limits(y = 0) + 
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0, 0.01)) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, c(0,NA))) + 
    ggplot2::labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
    ) + 
    ggplot2::theme_minimal()
  
  if (draw_boots) {
    p <- p + ggplot2::geom_line(data = dat$boot_wts, ggplot2::aes(x = p, y = wt, group = rep), color = boot_color, alpha = 100 * boot_alpha / R)
  }
  
  p + ggplot2::geom_line(aes(x = p, y = wt), color = color, linewidth = linewidth)
    
}
