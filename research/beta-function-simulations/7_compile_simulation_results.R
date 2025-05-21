library(tidyverse)
library(tictoc)
library(simhelpers)

params <- readRDS("simulation-beta-function/simulation_parameters.rds")

res_list <- tibble(
  file = list.files("simulation-beta-function/batch-results", pattern = "simulation_results_batch", full.names = TRUE)
) %>%
  mutate(
    row = str_extract(file, "batch[0-9]+.rds") |> str_sub(6,-5) |> as.integer()
  )
nrow(res_list)

outstanding_conditions <-
  params %>%
  anti_join(res_list, by = "row")

outstanding_conditions 

outstanding_conditions %>%
  select(row) %>%
  write_csv(file = "simulation-beta-function/batches-to-run.csv", col_names = FALSE)

tic()
sim_res_summary <- map_dfr(res_list$file, .f = readRDS)
sim_res_summary$file <- res_list$file
toc()

# file too big
write_rds(sim_res_summary, file = "simulation-beta-function/results/sim-beta-function-results.rds")
