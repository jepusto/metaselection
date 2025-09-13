library(metaselection)
library(metadat)
library(metafor)
library(tidyverse)
library(here)
#devtools::load_all()

setwd("C:/GitHub/metaselection/research/")

sci <- read_csv(here('tests/testdata','ScienceMeta.csv'))

sci <- sci |>
  rename(esid = es_id, studyid = study_id) |>
  mutate(
    d = yi,
    Var_d = vi,
    sd_d = sqrt(Var_d),
    var_d = Var_d,
    p_onesided = 1 - pnorm(d / sd_d)
  )

# Beta model - mean only
sci_adj_beta_mean <- selection_model(
  dat = sci,
  yi = d,
  sei = sd_d,
  pi = p_onesided,
  cluster = studyid,
  selection_type = "beta",
  steps = c(.025,.975),
  make_sandwich = TRUE
)

pvals <- seq(0.001, 0.999, .001)
selection_curves_beta_mean <- selection_wts(sci_adj_beta_mean, pvals = pvals)

sci_beta_mean <- ggplot(selection_curves_beta_mean) +
  aes(x = p, y = wt) +
  geom_vline(xintercept = c(.025, .975), linetype = "dashed") +
  geom_area(alpha = 0.5) +
  scale_x_continuous(
    limits = c(0,1),
    expand = expansion(0,0),
    transform = "asn",
    breaks = c(0,.025, .10, .25, .5, .75, .90, .975, 1)
  ) +
  theme_light() +
  labs(x = "p-value (one-sided)", y = "Selection weight") +
  theme(
    legend.position = "none",
    strip.text = element_text(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

ggsave("SREE2025-presentation/sci_beta_mean.png", plot = sci_beta_mean, width = 10, height = 6, dpi = 300)

# Beta model - moderator
sci_adj_beta_mod <- selection_model(
  dat = sci,
  yi = d,
  sei = sd_d,
  pi = p_onesided,
  cluster = studyid,
  selection_type = "beta",
  steps = c(.025,.975),
  make_sandwich = TRUE,
  mean_mods = ~ outcome_type_author
)

pvals <- seq(0.001, 0.999, .001)
selection_curves_beta_mod <- selection_wts(sci_adj_beta_mod, pvals = pvals)

sci_beta_mod <- ggplot(selection_curves_beta_mod) +
  aes(x = p, y = wt) +
  geom_vline(xintercept = c(.025, .975), linetype = "dashed") +
  geom_area(alpha = 0.5) +
  scale_x_continuous(
    limits = c(0,1),
    expand = expansion(0,0),
    transform = "asn",
    breaks = c(0,.025, .10, .25, .5, .75, .90, .975, 1)
  ) +
  theme_light() +
  labs(x = "p-value (one-sided)", y = "Selection weight") +
  theme(
    legend.position = "none",
    strip.text = element_text(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

ggsave("SREE2025-presentation/sci_beta_mod.png", plot = sci_beta_mod, width = 10, height = 6, dpi = 300)
