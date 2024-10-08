#' Interleaved Learning Meta-Analysis
#'
#' Meta-analytic dataset containing results of primary studies examining the 
#' effect of interleaved learning. 
#'
#' @format A data frame with 238 rows and 21 variables:
#' \describe{
#'   \item{study}{name of the study}
#'   \item{esid}{identifier for the effect size}
#'   \item{g}{effect size in form of Hedges' g}
#'   \item{vg}{corresponding variance of the effect size}
#'   \item{sample_id}{identifier for each sample within a study}
#'   \item{article_num}{identifier for the publication}
#'   \item{item_num}{code for items with 
#'    (1) Paintings including mostly impressionistic paintings of different artists coded as 0, 
#'    (2) Naturalistic photographs such as pictures of birds, butterflies coded as 3, 
#'    (3) Artificial pictures, that is, pictures of artificial objects or creatures coded as 1, 
#'    (4) Mathematical tasks, such as calculating the volume of geometric solids or the use of significance tests coded as 2, 
#'    (5) Expository Texts, which included plain expository texts and combinations of texts and other media in a multimedia/interactive media environment coded as 5, 
#'    (6) Words, such as names that belonged to different conceptual categories, pronunciation rules, or translations in different languages coded as 4, 
#'    (7) Tastes such as liquids with different tastes coded as 6}
#'   \item{age}{mean age of the participants}
#'   \item{grey_lit}{indicator for grey literature with 1 indicating theses and dissertations and 0 indicating articles and conference papers}
#'   \item{design}{indicator for research design with 1 indicating within-participants design and 0 indicating between-participants design}
#'   \item{students}{indicator type of sample with 1 indicating samples with only students and 0 indicating all other types of samples}
#'   \item{retention_interval}{indicator for retention interval length with 1 indicating long (>= 20 min) intervals and 0 indicating short (< 20 min) intervals}
#'   \item{intentionality}{indicator for intentional learning designs with 1 indicating incidental learning designs and 0 indicating intentional learning designs}
#'   \item{transfer_retention}{indicator for type of tests with 1 indicating transfer tests and 0 indicating retention tests}
#'   \item{simultaneity}{indicator for type of presentation with 1 indicating simultaneous presentation and 0 indicating successive presentation}
#'   \item{zmean_within}{standardized ratings of similarity within categories}
#'   \item{zmean_between}{standardized ratings of similarity between categories}
#'   \item{zmean_complex}{standardized ratings of complexity}
#'   \item{zmean_fam}{standardized ratings of familiarity}
#'   \item{zmean_cur}{standardized ratings of curiosity}
#'   \item{spaced}{indicator for spacing with spaced indicating designs with spaced items (10 or 30 sec time intervals), 
#'    nonspaced indicating designs with immediate succession of items (less than 2 sec), 
#'    -1 indicating designs with distractors between presentation of items are coded as distract and items not included in the analysis of spacing between items}
#' }
#'
#'@source \href{https://osf.io/7u253/}{OSF page for the project}
#'
#' @references Brunmair, M., & Richter, T. (2019). Similarity matters:
#'   A meta-analysis of interleaved learning and its moderators.  
#'   \emph{Psychological Bulletin, 145}(11), 1029.
#'   \doi{doi:10.1037/bul0000209}
#'
"interleaved_learning"