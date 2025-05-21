library(tidyverse)
library(tictoc)
library(simhelpers)

params <- readRDS("simulation-step-function/simulation_parameters.rds")

tic()
sim_res_summary <- 
  list.files("/Users/mjoshi/Library/CloudStorage/OneDrive-SharedLibraries-AIR/IES Selective Reporting - Documents/5_Simulation_Results", pattern = "simulation_results_batch", full.names = TRUE) %>%
  map_dfr(.f = readRDS)
toc()


params %>% 
  mutate(batchID = row_number()) %>%
  anti_join(sim_res_summary) %>%
  select(mean_smd:batch, batchID)

tic()
sim_res <-
  sim_res_summary %>%
  #slice(1:10) %>%
  select(-seed, -run_date) %>%
  unnest(res) %>% 
  ungroup() 
toc()




head(sim_res)
table(sim_res$param)
table(sim_res$estimator)
table(sim_res$iterations)
table(sim_res$R_conv)

sim_res %>%
  filter(str_detect(estimator, "ML|hybrid")) %>%
  group_by(estimator, R_conv) %>% 
  count()

sim_res %>%
  filter(is.na(R_conv)) %>%
  group_by(estimator) %>%
  count()


sim_res %>%
  filter(!str_detect(estimator, "ML|hybrid")) %>%
  group_by(estimator, R_conv) %>%
  count() %>%
  View()


sim_res_analyze <- 
  sim_res %>%
  mutate(true_smd = mean_smd) %>%
  filter(param == "beta") %>%
  filter(estimator == "ML" & R_conv == 0 | 
         str_detect(estimator, "hybrid") & R_conv <= 2 |
         !str_detect(estimator, "ML|hybrid"))

sim_res_analyze %>%
  group_by(estimator, R_conv) %>%
  count() %>%
  View()

# which R conv to take out?
perf_results_abs <- 
  sim_res_analyze %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, omega, steps, estimator) %>%
  group_modify(~ calc_absolute(.x, estimates = Est, 
                               true_param = true_smd,
                               perfm_criteria = c("bias", "rmse"))) %>%
  ungroup()

table(perf_results_abs$estimator)

perf_results_ci <- 
  sim_res_analyze %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, omega, steps, estimator) %>%
  group_modify(~ calc_coverage(.x, lower_bound = CI_lo,
                               upper_bound = CI_hi,
                               true_param = true_smd)) %>%
  ungroup()

table(perf_results_ci$estimator)

rate_conv <- 
  sim_res_analyze %>%
  filter(param == "beta") %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, omega, steps, estimator, R_conv) %>%
  summarize(n = n()) %>%
  mutate(p = n / sum(n)) 

# just exploring

table(perf_results_abs$m)



perf_results <- left_join(perf_results_abs, perf_results_ci)

table(perf_results$estimator)



# bias rmse and coverage --------------------------------------------------

# Bias 
graph_res <- perf_results %>% 
  mutate(weights = as.character(weights),
         m = as.character(m),
         weights = factor(weights, levels = c("0.05", "0.1", "0.2", "0.5", "1")),
         m = factor(m, levels = c("15", "30", "60", "90", "120", "200"))) %>%
  filter(!(estimator %in% c("PET", "PEESE"))) 

table(graph_res$estimator)



create_plot <- function(dat = graph_res, 
                        x_axis = weights,
                        y_axis,
                        y_int = 0,
                        low_coord = 0,
                        x_lab,
                        y_lab) {
  
  
  dat %>%
    ggplot(aes(x = {{x_axis}}, y = {{y_axis}}, color = estimator, fill = estimator)) +
    geom_boxplot(alpha = .5) +
    geom_hline(yintercept = y_int, linetype = "dashed") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    coord_cartesian(ylim = c(low_coord, 1)) +
    scale_x_discrete(labels = function(x) 
      stringr::str_wrap(x, width = 3))+
    facet_grid(tau ~ mean_smd, 
               labeller = label_bquote(rows = tau == .(tau),
                                       cols = mu == .(mean_smd)),
               scales = "free_y") +
    labs(x = x_lab, y = y_lab, fill = "Method", color = "Method") + 
    theme_bw() +
    theme(legend.position = "bottom")
  
  
}

create_plot(dat = graph_res,
            x_axis = weights,
            y_axis = bias,
            x_lab = "Selection Weights",
            y_lab = "Bias")

ggsave("simulation-step-function/results_preliminary/bias.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res %>% filter(str_detect(estimator, "ML|hybrid")),
            x_axis = weights,
            y_axis = bias,
            x_lab = "Selection Weights",
            y_lab = "Bias")

ggsave("simulation-step-function/results_preliminary/bias_rve.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res,
            x_axis = m,
            y_axis = bias,
            x_lab = "Number of Studies",
            y_lab = "Bias")

ggsave("simulation-step-function/results_preliminary/bias_num_studies.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res %>% filter(str_detect(estimator, "ML|hybrid")),
            x_axis = m,
            y_axis = bias,
            x_lab = "Number of Studies",
            y_lab = "Bias")

ggsave("simulation-step-function/results_preliminary/bias_rve_num_studies.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res,
            x_axis = weights,
            y_axis = rmse,
            x_lab = "Selection Weights",
            y_lab = "RMSE")

ggsave("simulation-step-function/results_preliminary/rmse.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res %>% filter(str_detect(estimator, "ML|hybrid")),
            x_axis = weights,
            y_axis = rmse,
            x_lab = "Selection Weights",
            y_lab = "RMSE")



ggsave("simulation-step-function/results_preliminary/rmse_rve.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res %>% filter(str_detect(estimator, "hybrid")),
            x_axis = weights,
            y_axis = rmse,
            x_lab = "Selection Weights",
            y_lab = "RMSE")


ggsave("simulation-step-function/results_preliminary/rmse_hybrid.png", device = "png", width = 12, height = 8, dpi = 500)

create_plot(dat = graph_res,
            x_axis = m,
            y_axis = rmse,
            x_lab = "Number of Studies",
            y_lab = "Bias")

ggsave("simulation-step-function/results_preliminary/rmse_num_studies.png", device = "png", width = 12, height = 8, dpi = 500)


# check what is going on with RMSE
check_rmse <- sim_res_analyze %>%
  filter(estimator == "ML") %>%
  mutate(check = (Est - mean_smd) ^ 2)

create_plot(dat = graph_res,
            x_axis = weights,
            y_axis = coverage,
            y_int = .95,
            low_coord = 0.6,
            x_lab = "Selection Weights",
            y_lab = "Coverage")

ggsave("simulation-step-function/results_preliminary/coverage.png", device = "png", width = 12, height = 8, dpi = 500)


create_plot(dat = graph_res %>% filter(str_detect(estimator, "ML|hybrid")),
            x_axis = weights,
            y_axis = coverage,
            y_int = .95,
            low_coord = 0.6,
            x_lab = "Selection Weights",
            y_lab = "Coverage")

ggsave("simulation-step-function/results_preliminary/coverage_rve.png", device = "png", width = 12, height = 8, dpi = 500)




create_plot(dat = graph_res,
            x_axis = m,
            y_axis = coverage,
            y_int = .95,
            low_coord = 0.6,
            x_lab = "Number of Studies",
            y_lab = "Coverage")

ggsave("simulation-step-function/results_preliminary/coverage_num_studies.png", device = "png", width = 12, height = 8, dpi = 500)


# hybrid ------------------------------------------------------------------

hybrid_dat <- sim_res %>%
  filter(str_detect(estimator, "hybrid|ML"), param == "beta") %>%
  filter(!(estimator == "ML" & R_conv == 1)) %>%
  filter(estimator != "hybrid-B")

glimpse(hybrid_dat)

hybrid_dat %>%
  group_by(rep, estimator, R_conv) %>% 
  count() %>%
  View()


hybrid_dat <- 
  hybrid_dat %>%
  mutate(order = case_when(estimator == "hybrid-A" & R_conv %in% c(1, 2) ~ 1,
                           estimator == "hybrid-C" & R_conv %in% c(1, 2) ~ 2, 
                           estimator == "hybrid-D" & R_conv %in% c(1, 2) ~ 3,
                           estimator == "ML" ~ 4,
                           TRUE ~ as.numeric(NA)))


hybrid_dat %>%
  select(rep, estimator, R_conv, Est, order) %>%
  View()

hybrid_dat_selected <- 
  hybrid_dat %>%
  filter(!is.na(order)) %>%
  group_by(mean_smd, tau, cor_mu, cor_sd, weights, m, omega, steps, rep) %>%
  filter(order == min(order)) %>%
  ungroup() 


hybrid_dat_selected %>%
  select(rep, estimator, R_conv, Est, order) %>%
  View()


check <- sim_res_analyze %>%
  filter(estimator == "ML")


check_miss <- anti_join(hybrid_dat_selected, check, by = c("mean_smd",
                                                           "tau",
                                                           "cor_mu",
                                                           "cor_sd",
                                                           "weights",
                                                           "m",
                                                           "iterations",
                                                           "omega",
                                                           "steps",
                                                           "rep"))

hybrid_dat_selected <- 
  hybrid_dat_selected %>%
  mutate(specific = estimator, 
         estimator = "hybrid")

glimpse(hybrid_dat_selected)

sim_res_clean <- 
  sim_res_analyze %>%
  filter(!str_detect(estimator, "hybrid")) %>%
  bind_rows(hybrid_dat_selected)


table(sim_res_clean$estimator)

glimpse(sim_res_clean)

# file too big
#write_rds(sim_res_clean, file = "simulation-step-function/results_preliminary/sim_res_clean.rds")
