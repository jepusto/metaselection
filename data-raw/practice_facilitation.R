library(tidyverse)

practice_facilitation_raw <- read_csv("data-raw/Baskerville et al. (2012).csv")

practice_facilitation <-
  practice_facilitation_raw %>%
  mutate(
    design = str_extract(`Trial Characteristics`, "Design: [RCT-]+") |>
      str_sub(9,-1),
    allocation_concealed = str_extract(`Trial Characteristics`, "Allocation concealed: [YN]") |>
      str_sub(23,-1),
    blinded = str_extract(`Trial Characteristics`, "Blindedb: [YN]") |>
      str_sub(11,-1),
    intent_to_treat = str_extract(`Trial Characteristics`, "Intent to treat: [YN]") |>
      str_sub(18,-1),
    across(
      c(allocation_concealed, blinded, intent_to_treat), 
      ~ recode(.x, "Y" = "Yes", "N" = "No")
    ),
    follow_up = str_extract(`Months Follow-up % retention`, "^[0-9]+") |> as.integer(),
    retention_pct = str_extract(`Months Follow-up % retention`, "\\([0-9.]+\\)") |> 
      str_sub(2,-2) |> 
      as.numeric(),
    SMD = str_extract(`Effect Size`, "^[0-9.]+"),
    SE = str_extract(`Effect Size`, "\\([0-9.]+\\)") |>
      str_sub(2, -2) |> 
      as.numeric()
  ) %>%
  select(
    author = Author, score = Score, 
    design, allocation_concealed, blinded, intent_to_treat,
    outcome = `Outcome Measure`, 
    follow_up, retention_pct, 
    SMD, SE
  )

practice_facilitation %>%
  View()


usethis::use_data(practice_facilitation, overwrite = TRUE)
