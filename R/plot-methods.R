#-------------------------------------------------------------------------------
# Calculate selection weights for specified p-values

#' @export

selection_wts <- function(mod, pts, ...) UseMethod("selection_wts")

#' @export

selection_wts.default <- function(mod, pts, params, steps) {
  mod <- match.arg(mod, c("step","beta"))
  if (mod == "step") {
    if (length(params) != length(steps)) stop("The length of `params` must be equal to the length of `steps`.")
    sel_fun <- step_fun(cut_vals = steps, weights = exp(params))
  } else if (mod == "beta") {
    if (length(params) != 2L) stop("The `params` argument must be a vector of length 2.")
    if (length(steps) != 2L) stop("The `steps` argument must be a vector of length 2.")
    sel_fun <- beta_fun(
      delta_1 = exp(params[1]), delta_2 = exp(params[2]),
      trunc_1 = steps[1], trunc_2 = steps[2]
    )
  }
  
  sel_fun(pts)
}

#' @export

selection_wts.step.selmodel <- function(mod, pts, bootstraps = TRUE) {
  
  if (!is.null(res$cl$sel_mods)) stop("selection_wts() is not available for models that include moderators of the selection parameters.")
  
  zeta <- mod$est$Est[grepl("zeta", mod$est$param)]
  
  sel_fun <- step_fun(cut_vals = mod$steps, weights = exp(zeta))
  
  wts <- data.frame(
    p = pts, 
    wt = sel_fun(pts)
  )
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)

  R <- mod$cl$R
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
    trunc_1 = mod$steps[1], trunc_2 = mod$steps[2]
  )
  
  wts <- data.frame(
    p = pts, 
    wt = sel_fun(pts)
  )
  
  if (!inherits(mod, "boot.selmodel") | !bootstraps) return(wts)
  
  R <- mod$cl$R
  sel_param_boots <- mod$bootstrap_reps[grepl("zeta", mod$bootstrap_reps$param),]
  sel_param_boots$rep <- rep(1:R, each = nrow(sel_param_boots) / R)
  sel_param_boots <- reshape(sel_param_boots, direction = "wide", idvar = "rep", timevar = "param")
  sel_param_boots$rep <- NULL
  
  boot_wts <- apply(
    sel_param_boots, 1, 
    \(zeta) beta_fun(
      delta_1 = exp(zeta[1]), delta_2 = exp(zeta[2]), 
      trunc_1 = mod$steps[1], trunc_2 = mod$steps[2]
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

select_plot <- function(mod, pts, ...) UseMethod("selection_plot")

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
  col = "blue",
  alpha = 0.5,
  ...
) {
  
  require(ggplot2)
  
  if (!is.null(res$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)

  dat <- selection_wts(mod, pts = pts, bootstrap = FALSE)
  
  ggplot(dat) + 
    aes(x = p, y = wt) + 
    expand_limits(y = 0) + 
    scale_x_continuous(limits = c(0, 1), expand = expansion(0, 0.01)) + 
    scale_y_continuous(expand = expansion(0, c(0,NA))) + 
    geom_area(fill = col, alpha = alpha) + 
    labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
    )
  
}

selection_plot.boot.selmodel <- function(
    mod, 
    pts = 200L, 
    col = "black",
    linewidth = 1.2, 
    boot_col = "blue",
    boot_alpha = 0.1,
    ...
) {
  
  require(ggplot2)
  
  if (!is.null(res$cl$sel_mods)) stop("selection_plot() is not available for models that include moderators of the selection parameters.")
  
  pts <- seq(0, 1, length.out = pts)
  
  dat <- selection_wts(mod, pts = pts, bootstraps = TRUE)
  R <- mod$cl$R
  
  ggplot(dat$wts) + 
    expand_limits(y = 0) + 
    geom_line(data = dat$boot_wts, aes(x = p, y = wt, group = rep), color = boot_col, alpha = 100 * boot_alpha / R) + 
    geom_line(aes(x = p, y = wt), color = col, linewidth = linewidth) + 
    scale_x_continuous(limits = c(0, 1), expand = expansion(0, 0.01)) + 
    scale_y_continuous(expand = expansion(0, c(0,NA))) + 
    labs(
      x = "p-value (one-sided)",
      y = "Selection weight"
    ) + 
    theme_minimal()
  
}
