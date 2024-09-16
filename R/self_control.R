#' What Works Clearinghouse sample size and effect size distribution data
#'
#' A dataset containing sample size and number of effect sizes seen in primary
#' studies evaluated by the What Works Clearinghouse.
#'
#' @format A tibble with 158 rows and 27 variables:
#' \describe{
#'   \item{studyid}{identifier for study}
#'   \item{esid}{identifier for effect sizes}
#'   \item{name}{name of the study}
#'   \item{g}{effect size in form of Hedges' g}
#'   \item{var_g}{variance of the effect size}
#'   \item{se_g}{standard error of the effect size}
#'   \item{outcome}{name of the outcome measure}
#'   \item{comparison}{conditions compared}
#'   \item{type_of_treatment}{type of treatment}
#'   \item{sample_population}{population examined in the study}
#'   \item{sample_age}{average age of the sample}
#' }
#'
#'
#'
#' @importFrom Rdpack reprompt
#'
"self_control"