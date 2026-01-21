# Data cleaning
library(tidyverse)
library(broom)
library(GGally)

# Parallel processing
library(future)
plan(multisession, workers = 8L)

devtools::load_all()

ed_dat <- readRDS("research/step-function-paper/ego-depletion/Carter-et-al-data.rds")

#-------------------------------------------------------------------------------
# Limit to published studies, exclude large outlier

ed_dat_clean <- 
  ed_dat %>%
  filter(
    pub_type == 1,
    yi < 2
  )

prior_spec <- define_priors(
  beta_mean = 0,
  beta_precision = 1 / 16, 
  beta_L = 4,
  tau_mode = 0.2, 
  tau_alpha = 1,
  lambda_mode = 0.5,
  lambda_precision = 1 / 54,
  lambda_L = 4
)

PML <- selection_model(
  data = ed_dat_clean, 
  yi = yi,
  vi = vi,
  cluster = study_id,
  selection_type = "step",
  steps = c(0.025, .5),
  estimator = "CML",
  priors = prior_spec,
  CI_type = "large-sample"
)

ARGL <- selection_model(
  data = ed_dat_clean, 
  yi = yi,
  vi = vi,
  cluster = study_id,
  selection_type = "step",
  steps = c(0.025, .5),
  estimator = "ARGL",
  priors = prior_spec,
  CI_type = "large-sample"
)

scores_PML <- step_score(theta = PML$est$Est, yi = ed_dat_clean$yi, sei = ed_dat_clean$sei, steps = c(.025, .500), priors = prior_spec, contributions = TRUE)
hessian_PML <- step_hessian(theta = PML$est$Est, yi = ed_dat_clean$yi, sei = ed_dat_clean$sei, steps = c(.025, .500), priors = prior_spec)
hess_inv_PML <- solve(hessian_PML)
score_j_PML <- rowsum(scores_PML %*% hess_inv_PML, group = ed_dat_clean$study_id)
crossprod(score_j_PML)

scores_ARGL <- step_hybrid_score(theta = ARGL$est$Est, yi = ed_dat_clean$yi, sei = ed_dat_clean$sei, steps = c(.025, .500), priors = prior_spec, contributions = TRUE)
hessian_ARGL <- step_hybrid_jacobian(theta = ARGL$est$Est, yi = ed_dat_clean$yi, sei = ed_dat_clean$sei, steps = c(.025, .500), priors = prior_spec)
hess_inv_ARGL <- solve(hessian_ARGL)
score_j_ARGL <- rowsum(scores_ARGL %*% t(hess_inv_ARGL), group = ed_dat_clean$study_id)
crossprod(score_j_ARGL)

neg_vals <- which(ed_dat_clean$yi < 0)
PML_dat <- as.data.frame(scores_PML) %>% mutate(neg_vals = ed_dat_clean$yi < 0)
ggpairs(PML_dat, mapping = aes(color = neg_vals))
ggpairs(score_j_ARGL)
ggpairs(scores_ARGL)


scores_ARGL[neg_vals,]
