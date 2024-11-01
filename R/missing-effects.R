#' @title Estimate number of missing effect sizes
#'
#' @description Estimate number of missing effect sizes based on a fitted \code{selmodel} object.
#'
#' @param object Fitted model of class inheriting from \code{"selmodel"}.
#'
#' @export
#'
#' @examples
#' 
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
#' missing_effects(res)
#'

missing_effects <- function(object) {
  if (!inherits(object, "selmodel")) stop("`object` must inherit from class 'selmodel'.")
  pr_obs <- object$predictions$Ai
  K_total <- sum(1 / pr_obs)
  K_missing <- K_total - object$n_effects
  data.frame(
    type = c("Observed","Missing","Total"),
    effects = c(object$n_effects, K_missing, K_total)
  )
}