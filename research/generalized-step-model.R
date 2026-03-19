library(tidyverse)

sel_probs <- 
  expand.grid(
    Aj = 0:10,
    psi = c(0, 0.5, 1, 2, 5, 10),
    lambda = c(.02, .05, .10, .20, .50, 1)
  ) %>%
  mutate(
    pr = 1 - (1 - lambda)^((Aj + 1)^psi),
    lambda_f = fct(as.character(lambda)) |> fct_rev()
  )

ggplot(sel_probs) + 
  aes(Aj, pr, color = lambda_f) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(0,10,2)) + 
  facet_wrap(~ psi, labeller = "label_both") + 
  theme_minimal() + 
  labs(
    x = expression(A[j]),
    y = expression(Pr(O[ij]==1 | A[j])),
    color = expression(lambda)
  )
