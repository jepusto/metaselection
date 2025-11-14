#' @title Calculate area under the selection weight function from a `selmodel`
#'   object
#'
#' @description Summarize the strength of selection by calculating the area
#'   under the selection weight function a `selmodel` object (excluding the area
#'   with weight fixed at 1). If the object has bootstrap replications, then a
#'   confidence interval will also be calculated.
#'
#' @param object Fitted model of class \code{"selmodel"}.
#' @param CI_type character string specifying the type of confidence interval to calculate, with options as in \code{"selection_model"}. If \code{NULL} (the default), it will be inherited from \code{object}.
#' @param conf_level desired coverage level for confidence intervals. If \code{NULL} (the default), it will be inherited from \code{object}, which has a default value of \code{.95}.
#' @param warn logical controlling whether warnings are displayed, with a default of \code{TRUE}.
#' 
#' @export
#'
#' @examples
#'
#' beta_noboot <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   selection_type = "beta",
#'   steps = c(0.025,0.5)
#' )
#'
#' p_area(beta_noboot)
#'
#' step_boot <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   selection_type = "step",
#'   steps = c(0.025,0.50),
#'   estimator = "ARGL",
#'   bootstrap = "multinomial",
#'   CI_type = "percentile",
#'   R = 6
#' )
#' 
#' p_area(step_boot)
#' 

p_area <- function(object, CI_type = NULL, conf_level = NULL, warn = TRUE) UseMethod("p_area")

#' @export

p_area.beta.selmodel <- function(object, CI_type = NULL, conf_level = NULL, warn = TRUE) {
  
  steps <- object$steps
  pdim <- object$param_dim
  zeta_index <- sum(pdim[1:2]) + 1:pdim[3]
  zeta <- object$est$Est[zeta_index]
  
  p_area <- calc_beta_area(zeta, steps)
  
  res <- data.frame(
    param = "p-area",
    Est = p_area
  )
  
  if (inherits(object, "boot.selmodel")) {
    CI_type <- get_p_area_CI_type(object, CI_type, warn = warn)
    if (is.null(conf_level)) conf_level <- object$conf_level
    
    area_CIs <- p_area_boot_CI(
      object, p_area = p_area, 
      CI_type = CI_type, conf_level = conf_level, 
      steps = steps, zeta_index = zeta_index
    )
    
    res <- cbind(res, area_CIs)
  }
  
  res
}

calc_beta_area <- function(zeta, steps) {
  lambda <- exp(zeta)
  pts <- stats::dbeta(steps, shape1 = lambda[1], shape2 = lambda[2])
  ints <- stats::pbeta(steps, shape1 = lambda[1], shape2 = lambda[2])
  (diff(ints) + (1 - steps[2]) * pts[2]) / (pts[1] * (1 - steps[1]))
}

#' @export

p_area.step.selmodel <- function(object, CI_type = NULL, conf_level = NULL, warn = TRUE) {
 
  if (!is.null(object$cl$sel_mods)) stop("p_area() is not available for models that include moderators of the selection parameters.")
  
  steps <- object$steps
  pdim <- object$param_dim
  zeta_index <- sum(pdim[1:2]) + 1:pdim[3]
  zeta <- object$est$Est[zeta_index]
  
  p_area <- calc_step_area(zeta, steps)
  
  res <- data.frame(
    param = "p-area",
    Est = p_area
  )
  
  if (inherits(object, "boot.selmodel")) {
    CI_type <- get_p_area_CI_type(object, CI_type, warn = warn)
    if (is.null(conf_level)) conf_level <- object$conf_level
    
    area_CIs <- p_area_boot_CI(
      object, p_area = p_area, 
      CI_type = CI_type, conf_level = conf_level, 
      steps = steps, zeta_index = zeta_index
    )
    
    res <- cbind(res, area_CIs)
  }
  
  res
  
}

calc_step_area <- function(zeta, steps) {
  lambda <- exp(zeta)
  widths <- diff(c(steps,1))
  sum(widths * lambda) / (1 - steps[1])
}

get_p_area_CI_type <- function(object, CI_type, warn = TRUE) {
  
  if (is.null(CI_type)) {
    CI_type <- eval(object$cl$CI_type)
  }
  
  if ("student" %in% CI_type) {
    CI_type <- setdiff(CI_type, "student")
    msg <- "`CI_type = 'studentized'` is not available for `p_area()`."
    if (length(CI_type) == 0L) {
      stop(msg)
    } else if (warn) {
      warning(msg)
    }
  }
  
  if ("BCa" %in% CI_type & is.null(object$jack_vals)) {
    CI_type <- setdiff(CI_type, "BCa")
    msg <- "`CI_type = 'BCa'` is only available when `object` includes `CI_type` of 'BCa'."
    if (length(CI_type) == 0L) {
      stop(msg)
    } else  if (warn) {
      warning(msg)
    }
  }
  
  CI_type
}

p_area_boot_CI <- function(object, p_area, CI_type, conf_level, steps, zeta_index) {

  f_area <- if (object$selection_type == "step") calc_step_area else calc_beta_area  
  
  area_boots <- 
    object$bootstrap_reps |>
    split(rep(1:object$R, each = sum(object$param_dim))) |>
    lapply(\(x) x$Est[zeta_index]) |>
    sapply(f_area, steps = steps)
  
  SE <- data.frame(SE = sd(area_boots))
  
  if ("BCa" %in% CI_type) {
    
    area_jack <- 
      object$jack_vals[,zeta_index,drop=FALSE] |> 
      apply(1, f_area, steps = steps)
    
    inf_vals <- p_area - area_jack
    
    CIs <- simhelpers::bootstrap_CIs(
      boot_est = area_boots, 
      est = p_area, 
      influence = inf_vals,
      CI_type = CI_type, level = conf_level
    )
    
  } else {
    
    CIs <- simhelpers::bootstrap_CIs(
      boot_est = area_boots, 
      est = p_area, 
      CI_type = CI_type, level = conf_level
    )
    
  }
  
  cbind(SE, CIs)
  
}