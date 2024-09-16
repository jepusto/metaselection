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
#'   \item{type_of_treatment}{type of treatment/training}
#'   \item{type_of_treatment}{length of treatment coded in days}
#'   \item{publication_status}{indicator for whether the study was published}
#'   \item{publication_status_new}{indicator for whether the study was published}
#'   \item{publication_year}{year that the study was published}
#'   \item{research_group}{indicator for whether researchers of the study belonged to the strength model group}
#'   \item{control_group_quality}{quality of the control group with active indicating that the control group worked on some task}
#'   \item{gender_ratio}{percentage of males in the sample}
#'   \item{type_of_outcome}{type of outcome}
#'   \item{subjectivity_of_outcome_measurement}{indicator for the subjectivity of the outcome measure}
#'   \item{lab_based_versus_real_world_behavior}{indicator for whether the behavior measured was assessed in lab or in the real-world}
#'   \item{stamina_versus_strength}{indicator for whether the outcome was assessed with (Stamina) or without a preceding effortful task (Strength)}
#'   \item{pre_test_measure}{indicator for thether there was a pre and post measure of outcome or only post measure}
#'   \item{sample_population}{population examined in the study}
#'   \item{sample_age}{average age of the sample}
#'   \item{attrition}{percentage of attrition}
#'   \item{partipant_compensation}{type of compensation received by the participants}
#'   \item{self_control_potential}{indicator for whether the outcome measure required utilization of full self-control potential}
#' }
#'
#'
#'
#' @importFrom Rdpack reprompt
#'
"self_control"