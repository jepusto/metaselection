#'Practice Facilitation Meta-Analysis
#'
#'Meta-analytic dataset containing results of primary studies examining the
#'effect of practice facilitation on the uptake of evidence-based practices
#'(EBPs) in primary care settings.
#'
#'@format A data frame with 23 rows and 3 variables:
#' \describe{
#'   \item{author}{First author and publication year of primary study report}
#'   \item{score}{Scored on a scale from 0 to 12, in which the higher the score, the higher the quality of the study methods.}
#'   \item{design}{Study design, with CCT = controlled clinical trial, C-RCT = cluster randomized controlled trial, RCT = randomized controlled trial.}
#'   \item{allocation_concealed}{Indicator for allocation concealment.}
#'   \item{blinded}{Indicator for whether study was single- or double-blinded.}
#'   \item{intent_to_treat}{Indicator for whether study adhered to intent-to-treat principle.}
#'   \item{outcome}{Description of outcome measure.}
#'   \item{follow_up}{Months of follow-up.}
#'   \item{retention_pct}{Percentage of sample retained at follow-up.}
#'   \item{SMD}{effect size in form of Hedges' g}
#'   \item{SE}{corresponding variance of the effect size}
#' }
#'
#'@source \href{https://doi.org/10.1370/afm.1312}{Table 1 of Baskerville et al.
#'  (2012)}
#'
#'@references Baskerville, N. B., Liddy, C., & Hogg, W. (2012). Systematic
#'  review and meta-analysis of practice facilitation within primary care
#'  settings. \emph{Annals of Family Medicine, 10}(1), 63-74.
#'  \doi{doi:10.1370/afm.1312}
#'
"practice_facilitation"