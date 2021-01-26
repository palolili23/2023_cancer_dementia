library(rio)
library(lubridate)
library(tidyverse)
library(here)
library(magrittr)
# library(summarytools)

# Import first cohort -----------------------------------------------------

data <- import(here::here("01_data", "no_ltfu", "wide_noltfu.RData"))

data_long <- data %>% 
  select(id, cohort, dementia, cancer, cancer_prev, cancer_date, dementia_date, death_2015, starts_with("t2"),
         end_fup_2015, e1,
         sex, age_0, education_three_levels, apoe4, hd_prev, hd_date, diabetes_prev, diabetes_date, stroke_prev, stroke_date) %>% 
  group_by(id) %>% 
  slice(rep(1:n(), each = t2death_y + 1)) %>% #repeat rows by #years until death
  ungroup() %>%
  # select(everything(), - matches("[1,2,3,4,5]$"), apoe4) %>% ## delete all covariates
  group_by(id) %>%
  mutate(year = year(e1),
         year = year + (row_number() - 1)) %>% #an easy roll for years + 1
  ungroup()

### Create the repeated measurement data separate, by visit
rep_cov <- data %>% 
  select(id, 
         starts_with("e"),
         starts_with("bmi"), 
         starts_with("ht"), 
         starts_with("chol"),
         starts_with("oh"),
         starts_with("db_med"),
         starts_with("smoke"),
         starts_with("ht_drug"),
         starts_with("sbp"),
         starts_with("hdl"),
         starts_with("cesd"),
         -c(end_fup_2015, education_three_levels, cesd6, evt_reg)) %>% 
  pivot_longer(-1,
               names_to = c(".value", "visit"),
               names_sep = -1,
               values_drop_na = TRUE) %>%
  arrange(id, visit) %>% 
  drop_na(e) %>% #some have values on covariates but no date of visit
  mutate(year = year(e))

### check real missing data on covariates before carrying forward
rep_cov %>% 
  select(-id, - year, - e) %>% 
  summarise_all(~sum(is.na(.x))) 

rep_cov %>% 
  select(-id, - year, - e) %>% 
  summarise_all(~sum(is.na(.x))/n()) 

#### merge and fill
var_names <- colnames(rep_cov)

data_long <- data_long %>%
  left_join(rep_cov, by = c("id", "year")) %>%
  group_by(id) %>%
  fill(var_names[-1]) %>% 
  ungroup()

# data_long %>% 
#   VIM::aggr(plot = FALSE)
# 

# Generate visit and comorbidities variables ------------------------------------------------

library(zoo)

data_long <- data_long %>% 
  group_by(id) %>% 
  mutate(v1_bin = 1,
         v2_bin = ifelse((visit == 2), 1, NA),
         v2_bin = na.locf(v2_bin, na.rm =  FALSE),
         v2_bin = ifelse(is.na(v2_bin), 0, v2_bin),
         v3_bin = ifelse((visit == 3), 1, NA),
         v3_bin = na.locf(v3_bin, na.rm = FALSE),
         v3_bin = ifelse(is.na(v3_bin), 0, v3_bin),
         v4_bin = ifelse((visit == 4), 1, NA),
         v4_bin = na.locf(v4_bin, na.rm = FALSE),
         v4_bin = ifelse(is.na(v4_bin), 0, v4_bin)) %>% 
  ungroup()


data_long <- data_long %>% 
  group_by(id) %>% 
  mutate(
    hd_v = ifelse(hd_prev == 1, 1, NA),
    hd_date = as.numeric(year(hd_date)), 
    hd_v = ifelse((!is.na(hd_date) & hd_date == year), 1, hd_v),
    hd_v = na.locf(hd_v, na.rm = FALSE),
    hd_v = ifelse(is.na(hd_v), 0, hd_v),
    diab_v = ifelse(diabetes_prev == 1, 1, NA),
    diabetes_date = as.numeric(year(diabetes_date)), 
    diab_v = ifelse((!is.na(diabetes_date) & diabetes_date == year), 1, diab_v),
    diab_v = na.locf(diab_v, na.rm = FALSE),
    diab_v = ifelse(is.na(diab_v), 0, diab_v),
    stroke_v = ifelse(stroke_prev == 1, 1, NA),
    stroke_date = as.numeric(year(stroke_date)), 
    stroke_v = ifelse((!is.na(stroke_date) & stroke_date == year), 1, stroke_v),
    stroke_v = na.locf(stroke_v, na.rm = FALSE),
    stroke_v = ifelse(is.na(stroke_v), 0, stroke_v),
  )%>% 
  ungroup()

data_long %<>% 
  group_by(id) %>% 
  mutate(
    cancer_v = ifelse(cancer_prev == 1, 1, NA),
    cancer_date = as.numeric(year(cancer_date)),
    cancer_v = ifelse((!is.na(cancer_date) & cancer_date == year), 1, cancer_v),
    cancer_v = na.locf(cancer_v, na.rm = FALSE),
    cancer_v = ifelse(is.na(cancer_v),0, cancer_v)
  )

cancer_patients <- data %>% 
  filter(cancer == 1) %>% pull(id)

data_long %>%
  filter(id %in% cancer_patients) %>%
  select(cancer, cancer_v, cancer_date, year, dementia, death_2015) %>% View()
### Some checks

id_years <- data_long %>% 
  group_by(id) %>% 
  tally() %>% select(id)

test <- data %>% anti_join(id_years, by = "id")

#### until know everyone is followed until t2death
### To plot, I want to keep those who had dementia until they died

### Another dataset should be until dementia time

data_long %<>%
  group_by(id) %>% 
  arrange(id, year) %>% 
  mutate(time = row_number(),
         exposure = ifelse(cancer == 0, 0, NA),
         exposure = ifelse((t2cancer_y == time | t2cancer_y == 0) & cancer == 1, 1, exposure),
         exposure = zoo::na.locf(exposure, na.rm = FALSE),
         exposure = ifelse(is.na(exposure), 0, exposure)) %>% 
  ungroup() %>% 
  select(id, year, max, time, exposure, dementia, everything())


####KIM: NEXT LINES DELETES THOSE WITH T2DEATH_2015_Y==0 (this variable is already set at 0 in syntax 02_clean_long_kim_adjusted)
#### but max needs to be set at 1 here instead of 0

data_long$max <- ifelse(data_long$max==0&data_long$t2death_2015_y==1,1,data_long$max)

data_long %<>%
  group_by(id) %>%
  mutate(
    outcome_plr = ifelse(time == max, outcome, 0),
    competing_plr = ifelse(time == max, competing, 0),
    no_cr = 1 - competing_plr,
    cens_plr = ifelse(time == max, cens, 0),
    no_cens = 1 - cens_plr,
  ) %>%
  filter(time <= max) %>%
  ungroup() %>%
  arrange(id, time)

data_long_no_cens <- data_long %>% 
  mutate(
    outcome_plr = case_when(
      cens_plr == 1 ~ NA_real_,
      competing_plr == 1 ~ NA_real_,
      TRUE ~ outcome_plr),
    competing_plr = ifelse(cens_plr == 1, NA, competing_plr),
    no_cr = ifelse(is.na(competing_plr), NA, no_cr))


export(data_long, here::here("01_data","no_ltfu", "long_noltfu.RData"))
