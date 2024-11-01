#' @title Predict method for fitted `step.selmodel`
#'
#' @description Calculate predicted values for the mean (location),
#'   heterogeneity (scale), and selection parameters from a fitted model of
#'   class \code{"step.selmodel"}.
#'
#' @param object Fitted model of class inheriting from \code{"step.selmodel"}.
#' @param newdata Optionally, a data frame containing observations for which to
#'   generate predictions. If \code{NULL} (the default), predictions are
#'   generated for the data frame used to fit the model.
#' @param location logical with `TRUE` (the default) indicating that the
#'   predicted value of the mean effect size should be returned.
#' @param scale logical with `TRUE` (the default) indicating that the predicted
#'   value of the heterogeneity should be returned.
#' @param selection logical with `TRUE` (the default) indicating that the
#'   predicted selection parameters should be returned.
#' @param obs_prob logical with `TRUE` (the default) indicating that the
#'   predicted probability of observing an effect size should be returned.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
#' @examples
#' res <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   selection_type = "step",
#'   steps = 0.025,
#'   mean_mods = ~ 0 + sample_population,
#'   var_mods = ~ 0 + sample_population
#' )
#'
#' predict(res)
#'
#' newdat <- data.frame(sample_population = c("Students","Students"))
#' predict(res, newdata = newdat)



predict.step.selmodel <- function(
  object, 
  newdata = NULL, 
  location = TRUE,
  scale = TRUE,
  selection = TRUE, 
  obs_prob = TRUE,
  ...
) {
  
  if (is.null(newdata)) {
    
    predictions <- object$predictions
    
    res <- data.frame(row.names = row.names(object$mf))
    
  } else {

    cl <- object$cl
    sei <- eval(cl$sei, envir = newdata)
    steps <- object$steps
    
    X <- if (is.null(cl$mean_mods)) NULL else do.call(model.matrix, list(object = cl$mean_mods, data = newdata))
    U <- if (is.null(cl$var_mods)) NULL else do.call(model.matrix, list(object = cl$var_mods, data = newdata))
    Z0 <- if (is.null(cl$sel_zero_mods)) NULL else do.call(model.matrix, list(object = cl$sel_zero_mods, data = newdata))
    if (is.null(cl$sel_mods)) {
      Z <- NULL
    } else {
      if (is.list(cl$sel_mods)) {
        if (length(cl$sel_mods) != length(steps)) stop("sel_mods must be a list with length equal to the number of steps.")
        Z <- lapply(cl$sel_mods, model.matrix, data = newdata)
      } else {
        Z1 <- do.call(model.matrix, list(object = cl$sel_mods, data = newdata))
        Z <- rep(list(Z1), times = length(steps))
      }
    }
    
    predictions <- parse_step_params(
      theta = object$est$Est,
      sei = sei, pi = NULL, ai = NULL, 
      steps = steps, 
      X = X, U = U, Z0 = Z0, Z = Z, calc_Ai = TRUE
    )
    
    res <- data.frame(row.names = row.names(newdata))
  }
  
  
  if (location) res$mu <- predictions$mu
  if (scale) res$tau2 <- predictions$tausq
  if (selection) {
    colnames(predictions$lambda_full) <- paste0("lambda", 1:ncol(predictions$lambda_full) - 1L)
    res <- cbind(res, predictions$lambda_full)
  }
  if (obs_prob) res$obs_prob <- predictions$Ai
  
  return(res)
  
}

#' @title Predict method for fitted `beta.selmodel`
#'
#' @description Calculate predicted values for the mean (location) and
#'   heterogeneity (scale) from a fitted model of
#'   class \code{"beta.selmodel"}.
#'
#' @param object Fitted model of class inheriting from \code{"beta.selmodel"}.
#' @inheritParams predict.step.selmodel
#' 
#' @export
#'
#' @examples
#' res <- selection_model(
#'   data = self_control,
#'   yi = g,
#'   sei = se_g,
#'   cluster = studyid,
#'   selection_type = "beta",
#'   mean_mods = ~ 0 + sample_population,
#' )
#'
#' predict(res)
#'
#' newdat <- data.frame(sample_population = c("Students","Students"))
#' predict(res, newdata = newdat)

predict.beta.selmodel <- function(
    object, 
    newdata = NULL, 
    location = TRUE,
    scale = TRUE,
    obs_prob = TRUE,
    ...
) {
  
  if (is.null(newdata)) {
    
    predictions <- object$predictions
    res <- data.frame(row.names = row.names(object$mf))
    
  } else {
    
    cl <- object$cl
    sei <- eval(cl$sei, envir = newdata)
    steps <- object$steps
    
    X <- if (is.null(cl$mean_mods)) NULL else do.call(model.matrix, list(object = cl$mean_mods, data = newdata))
    U <- if (is.null(cl$var_mods)) NULL else do.call(model.matrix, list(object = cl$var_mods, data = newdata))

    predictions <- parse_beta_params(
      theta = object$est$Est,
      sei = sei, pi = NULL, 
      alpha = steps, 
      X = X, U = U, calc_Ai = TRUE
    )
    
    res <- data.frame(row.names = row.names(newdata))
    
  }
  
  if (location) res$mu <- predictions$mu
  if (scale) res$tau2 <- predictions$tausq
  if (obs_prob) res$obs_prob <- predictions$Ai
  
  return(res)
  
}