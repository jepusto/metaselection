library(tidyverse)
library(janitor)

#-------------------------------------------------------------------------------
# Read in data, combine files

setwd("research/step-function-paper")

load("ego-depletion/DataFrames/FC.Rdata")
load("ego-depletion/DataFrames/HG.Rdata")
load("ego-depletion/DataFrames/IA.Rdata")
load("ego-depletion/DataFrames/IP.Rdata")
load("ego-depletion/DataFrames/PA.Rdata")
load("ego-depletion/DataFrames/S.Rdata")
load("ego-depletion/DataFrames/ST.Rdata")
load("ego-depletion/DataFrames/WM.Rdata")
# 
# HG <- HG %>% 
#   mutate(
#     Exp = as.factor(as.character(Exp)),
#     Year = as.factor(as.character(Year))
#   )
# 
# IA <- IA %>% mutate(Year = as.factor(as.character(Year)))
# 
# PA <- PA %>% 
#   mutate(
#     Exp = as.factor(as.character(Exp)),
#     Year = as.factor(as.character(Year))
#   )
# 
# S <- S %>% mutate(Exp = as.factor(as.character(Exp)))
# 
# ST <- ST %>% mutate(Year = as.factor(as.character(Year)))
# 
# WM <- WM  %>% mutate(Year = as.factor(as.character(Year)))

ed_dat <- 
  bind_rows(
    `Food consumed` = FC,
    `Hand grip` = HG,
    `Impossible anagrams` = IA,
    `Impossible puzzles` = IP,
    `Possible anagrams` = PA,
    `Stroop` = S,
    `Standardized test` = ST,
    `Working memory` = WM,
    .id = "outcome"
  ) %>%
  clean_names() %>%
  mutate(
    study = paste(author, year),
    sample = paste(author, year, exp)
  ) %>%
  group_by(study) %>%
  mutate(study_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(es_id = row_number()) %>%
  select(
    study, study_id, sample, es_id, 
    outcome, author, exp, year, id, yi = g, vi = g_v,
    everything()
  )

ed_dat %>%
  mutate(outcome = "Combined") %>%
  bind_rows(ed_dat) %>%
  mutate(
    outcome = fct(outcome) |> fct_relevel("Combined")
  ) %>%
  group_by(outcome) %>%
  summarize(
    studies = n_distinct(study),
    samples = n_distinct(sample),
    effects = n(),
    studies_pub = n_distinct(study[pub_type == 1]),
    samples_pub = n_distinct(sample[pub_type == 1]),
    effects_pub = sum(pub_type == 1),
  )


saveRDS(ed_dat, file = "ego-depletion/Carter-et-al-data.rds")
