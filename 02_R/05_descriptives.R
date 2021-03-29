library(survival)
library(tidyverse)

data_long <-
  rio::import(here::here("01_data", "clean_data", "dementia_long.RData"))

data_wide <-
  data_long %>%
  group_by(id) %>% 
  slice(max(row_number())) %>% 
  ungroup()

## Total n
total_n <- dim(data_wide)[1]

## Counts

total_cancer <- data_wide %>% 
  count(cancer_v) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(cancer_count = paste0(prop, "%", " (n = ", n, ")")) %>% 
  select(cancer_count, cancer_v)
  
no_cancer <- total_cancer %>% 
  filter(cancer_v == 0) %>% 
  pull(1)

yes_cancer <- total_cancer %>% 
  filter(cancer_v == 1) %>% 
  pull(1)

total_status <- data_wide %>%
  group_by(cancer_v) %>% 
  count(status) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(total = paste0(prop, "%", " (n = ", n, ")")) %>% 
  select(cancer_v, status, total)

cancer_dem <- total_status %>% filter(cancer_v == 1, status == "Cancer and Dementia") %>% pull(3)
cancer_death <- total_status %>% filter(cancer_v == 1, status == "Dead") %>% pull(3)
cancer_alive <- total_status %>% filter(cancer_v == 1, status == "Only cancer") %>% pull(3)

no_cancer_dem <- total_status %>% filter(cancer_v == 0, status == "Only dementia") %>% pull(3)
no_cancer_death <- total_status %>% filter(cancer_v == 0, status == "Dead") %>% pull(3)
no_cancer_alive <- total_status %>% filter(cancer_v == 0, status == "Alive, free of cancer and dementia") %>% pull(3)

cancer_time <- data_long %>% 
  filter(cancer_v == 1 & lag(cancer_v) == 0) %>% 
  summarize(m = median(age_scale),
            l = quantile(age_scale, 0.25),
            u = quantile(age_scale, 0.75)) %>% 
  mutate(temp = paste0(m, " (IQR: ", l, "-", u, ")")) %>% 
  pull(temp)

dem_total <- data_wide %>% count(outcome_plr) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(temp = paste0(prop, "%", " (n = ", n, ")")) %>%
  filter(outcome_plr == 1) %>% 
  pull(temp)
  
dem_time <- data_long %>% 
  filter(dem_v == 1 & lag(dem_v) == 0) %>% 
  summarize(m = median(age_scale),
            l = quantile(age_scale, 0.25),
            u = quantile(age_scale, 0.75)) %>% 
  mutate(temp = paste0(m, " (IQR: ", l, "-", u, ")")) %>% 
  pull(temp)


## Baseline covs

mean_age <- data_wide %>% 
  summarize(a = round(mean(age_0),  2),
            b = round(sd(age_0),  2)) %>%
  mutate(temp = paste0(a, " (SD: ", b, ")")) %>% 
  pull(temp)


women_prop <- data_wide %>%
  count(sex) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(total = paste0(prop, "%", " (n = ", n, ")")) %>% 
  filter(sex == 1) %>% 
  pull(total)
  
# survminer::ggcompetingrisks(t2cancer_km) +
#   ggthemes::scale_fill_tableau() +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom")

## Probability of cancer and dementia at 5 years
data_wide <- data_wide %>% 
  mutate(cancer_death = ifelse(cancer_v == 0 & competing_plr == 1, 2, cancer_v))

t2cancer_km <-
  survfit(Surv(t2cancer_y, as.factor(cancer_death)) ~ 1, data_wide)


cancer_death_5y <- t2cancer_km %>% 
  broom::tidy() %>% 
  filter(state %in% c(1,2)) %>% 
  filter(between(time, 5.00, 5.01)) %>% 
  select(estimate, conf.low, conf.high, state) %>% 
  mutate_if(is.numeric, ~.*100) %>% 
  mutate_if(is.numeric, round, 1) %>%  
  mutate(est_print = paste0(estimate, " (", conf.low, ", ", conf.high, ")")) %>% 
  select(est_print, state)

death_5y <- cancer_death_5y %>% filter(state == 2) %>% pull(1)
cancer_5y <- cancer_death_5y %>% filter(state == 1) %>% pull(1)

