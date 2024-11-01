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
    
    res <- data.frame(row.names = row.names(object$mf))
    
    if (location) res$mu <- object$predictions$mu
    if (scale) res$tau2 <- object$predictions$tausq
    if (selection) {
      colnames(object$predictions$lambda_full) <- paste0("lambda", 1:ncol(object$predictions$lambda_full) - 1L)
      res <- cbind(res, object$predictions$lambda_full)
    }
    if (obs_prob) res$obs_prob <- object$predictions$Ai

    return(res)
    
  } else {
    
  }
  
  res
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
    
    res <- data.frame(row.names = row.names(object$mf))
    
    if (location) res$mu <- object$predictions$mu
    if (scale) res$tau2 <- object$predictions$tausq
    if (obs_prob) res$obs_prob <- object$predictions$Ai
    
    return(res)
    
  } else {
    
  }
  
  res
}