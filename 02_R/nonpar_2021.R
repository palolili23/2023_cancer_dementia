library(ggplot2)
library(tidyverse)
library(rio)
library(wesanderson)


# Import data -------------------------------------------------------------

data_wide <-
  import(here::here("01_data", "clean_data", "wide_noltfu.RData"))

data_long <-
  import(here::here("01_data", "clean_data", "dementia_long.RData"))

# 
data_wide %>% count(cancer, cancer_20, cancer_prev, cancer_inc)

## Use these numbers for chart
data_wide %>% count(cancer_20, dementia_20)

## idea of how numbers look like over follow-up
long_in_wide <- data_long %>% 
  group_by(id) %>% 
  slice(n()) %>% 
  ungroup()

long_in_wide %>% 
  count(cancer_v, outcome_plr, competing_plr,)

# id_a <- data_wide %>% filter(cancer_20 == 1) %>% pull(id)
# 
# id_b <- long_in_wide %>% filter(cancer_v == 1) %>% pull(id)
# 
# id_dif <- setdiff(id_b, id_a)

# no_cancer <- data_wide %>% 
#   filter(id %in% id_dif)
# 
# no_cancer %>% 
#   select(id, year, time, end_dementia_death, cancer_date, dementia_date,
#          end_fup_2015, cancer_v, outcome_plr, status) %>% 
#   view()


# Non parametric estimation, ever vs never, incident and ipwtd ------------


