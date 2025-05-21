
#-------------------------------------------------------------------------------
# Source packages and functions
library(dplyr)
library(tidyr)
library(purrr)
library(tictoc)
library(metafor)
library(clubSandwich)
library(metaselection)
progressr::handlers(global = TRUE)

source("research/beta-function-simulations/2_estimation.R")
source("research/beta-function-simulations/3_performance_criteria.R")
source("research/beta-function-simulations/4_simulation_driver.R")

tic()

check <-
  run_sim(
    iterations = 5L,
    mean_smd = 0.1,
    tau = 0.3,
    omega = 0,
    m = 60,
    cor_mu = 0.4,
    cor_sd = 0.05,
    delta_1 = 0.2,
    delta_2 = 0.9,
    step_models = c("3PSM","4PSM"),
    comparison_methods = "All",
    bootstrap = "two-stage",
    R_beta = c(49,99,199,299),
    summarize_performance = FALSE,
    seed = 20250541L
  )

toc()

true_params <- data.frame(
  param = c("beta", "gamma", "zeta1", "zeta2"),
  true_param = c(0.1, log(0.3^2), log(0.2), log(0.9))
)


check %>%
  merge(true_params) %>%
  calc_performance(winz = 2.5)
