library(metaselection)
library(ggplot2)
library(rddensity)

n_ES_fixed <- function(m) {
  data.frame(
    n = rep(60L, m),
    n_ES = rep(8L, m)
  )
}

dat <- r_meta(
  mean_smd = 0.2,
  tau = 0.1,
  omega = 0.05,
  m = 1e4,
  cor_mu = 0.4,
  cor_sd = 0.05,
  censor_fun = step_count_fun(cut_val = .025, weight = 0.2, psi = 0),  
  n_ES_sim = n_ES_fixed,
  m_multiplier = 5L
)


tstat_density <- rddensity(dat$t_i, c = qnorm(.975))
rdplotdensity(
  tstat_density, 
  X = dat$t_i,
  plotRange = range(dat$t_i)
)

dat %>%
  mutate(
    sig = if_else(p_onesided < .025, "A","N")
  ) %>%
ggplot() + 
  aes(d, color = sig, fill = sig) + 
  geom_density(trim = TRUE)
