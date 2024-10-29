#' @export

logLik.selmodel <- function (object, REML = FALSE, ...) {
  p <- sum(object$param_dim)
  N <- object$n_effects
  val <- object$ll
  attr(val, "nobs") <- N
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}
