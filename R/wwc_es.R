#' What Works Clearinghouse sample size and effect size distribution data
#'
#' A dataset containing sample size and number of effect sizes seen in primary
#' studies evaluated by the What Works Clearinghouse.
#'
#' @format A tibble with 615 rows and 2 variables:
#' \describe{
#'   \item{n}{the sample size of the primary study}
#'   \item{n_ES}{number of effect sizes in the primary study}
#' }
#'
#' @importFrom Rdpack reprompt
#'
#' @references 
#' \insertRef{wwc}{metaselection}
#'
"wwc_es"