library(rio)
library(lubridate)
library(tidyverse)
library(magrittr)


# Import data-----------------------------------------------------

data <- import(here::here("01_data", "clean_data", "wide_noltfu.RData"))

# Create long, as many rows as t2death in years + 1

data_long <- data %>% 
  select(id, cohort, dementia, cancer, cancer, cancer_prev, cancer_date,
         dementia_date, death_2015, starts_with("t2"),
         end_fup_2015, e1,
         sex, age_0, education, apoe4, 
         hd_prev, hd_date, diabetes_prev, diabetes_date, stroke_prev, 
         stroke_date) %>% 
  mutate(t2death_y = round(t2death_y,0)) %>% 
  filter(t2death_y >= 1) %>% 
  group_by(id) %>% 
  slice(rep(1:n(), each = t2death_y + 2)) %>% #repeat rows by #years until death
  mutate(year = year(e1),
         year = year + (row_number() - 1)) %>% #an easy roll for years + 1
  ungroup()


## Completes with rows for all missing years (everyone until 2016 because I extended an extra row)
data_long_2015 <- data_long %>% 
  complete(id, year)

### Make dataset with visit data, and transform to long
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
         -c(end_fup_2015, education, cesd6, evt_reg)) %>% 
  pivot_longer(-1,
               names_to = c(".value", "visit"),
               names_sep = -1,
               values_drop_na = TRUE) %>%
  arrange(id, visit) %>% 
  drop_na(e) %>% #some have values on covariates but no date of visit
  mutate(year = year(e))

## join and fill missign values (fills gaps between 2 visits by carrying forward)

data_long_2015 %<>% 
  left_join(rep_cov, by = c("id", "year"))

var_names <- colnames(data_long_2015)

data_long_2015 %<>% 
  group_by(id) %>%
  fill(var_names[-1]) %>% 
  ungroup()

### remove rows prior to study entry
data_long_2015 %<>% 
  group_by(id) %>% 
  filter(year >= year(e1)) %>% 
  ungroup()


# Make indicator of disease -----------------------------------------------
library(zoo)

data_long_2015 %<>% 
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


data_long_2015 %<>% 
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

data_long_2015 %<>% 
  group_by(id) %>% 
  mutate(
    cancer_v = ifelse(cancer_prev == 1, 1, NA),
    cancer_date = as.numeric(year(cancer_date)),
    cancer_v = ifelse((!is.na(cancer_date) & cancer_date == year), 1, cancer_v),
    cancer_v = na.locf(cancer_v, na.rm = FALSE),
    cancer_v = ifelse(is.na(cancer_v),0, cancer_v),
    dem_date = as.numeric(year(dementia_date)),
    dem_v = ifelse((!is.na(dementia_date) & dem_date == year), 1, NA),
    dem_v = na.locf(dem_v, na.rm = FALSE),
    dem_v = ifelse(is.na(dem_v),0, dem_v),
    death_year = as.numeric(year(end_fup_2015)),
    death_v = ifelse(death_year == year & death_2015 == 1, 1, NA),
    death_v = na.locf(death_v, na.rm = FALSE),
    death_v = ifelse(is.na(death_v),0, death_v)) %>% 
  ungroup()

data_long_2015 %<>% 
  group_by(id) %>% 
  mutate(cancer_lag = lag(cancer_v),
         cancer_lag = ifelse(is.na(cancer_lag), 0, cancer_lag),
         stroke_lag = lag(stroke_v),
         stroke_lag = ifelse(is.na(stroke_lag), 0, stroke_lag),
         hd_lag = lag(hd_v),
         hd_lag = ifelse(is.na(hd_lag), 0, hd_lag),
         diab_lag = lag(diab_v),
         diab_lag = ifelse(is.na(diab_lag), 0, diab_lag),
         sbp_lag = lag(sbp),
         sbp_lag = ifelse(is.na(sbp_lag), sbp, sbp_lag),
         bmi_lag = lag(bmi),
         bmi_lag = ifelse(is.na(bmi_lag), bmi, bmi_lag),
         smoke_lag = lag(smoke),
         smoke_lag = ifelse(is.na(smoke_lag), smoke, smoke_lag),
         ht_lag = lag(ht),
         ht_lag = ifelse(is.na(ht_lag), ht, ht_lag),
         ht_drug_lag = lag(ht_drug),
         ht_drug_lag = ifelse(is.na(ht_drug_lag), ht_drug, ht_drug_lag),
         ) %>% 
  ungroup()

data_long_2015 %>%
  select(
    id,
    cohort,
    year,
    e1,
    cancer,
    cancer_v,
    cancer_date,
    dementia,
    dem_v,
    dem_date,
    death_2015,
    death_v,
    death_year,
    t2death_y
  ) %>% View()

#### Create a "status" variable

data_long_2015 %<>% 
  mutate(status = case_when(
    cancer_v == 1 & dem_v == 0 ~ "Only cancer",
    cancer_v == 0 & dem_v == 1 ~ "Only dementia",
    cancer_v == 1 & dem_v == 1 ~ "Cancer and Dementia",
    TRUE ~ "Alive, free of cancer and dementia"
  ),
  status = ifelse(death_v == 1 & dem_v == 0, "Dead", status))

## Age scale

data_long_2015 %<>%
  group_by(id) %>% 
  mutate(age_scale = round(age_0,0) + (row_number() - 1)) %>% 
  ungroup()

## Export

export(data_long_2015, here::here("01_data", "clean_data", "data_long_2015.RData")) 


# Rows up to dementia diagnosis -------------------------------------------

data_long_dem <- data_long_2015 %>%
  group_by(id) %>% 
  arrange(id, year) %>% 
  mutate(end_dementia_death = ifelse(!is.na(dementia_date), 
                                     as.numeric(year(dementia_date)),
                                     as.numeric(year(end_fup_2015))),
         time = row_number() - 1) %>%
  filter(year <= end_dementia_death,
         time <= 20) %>%
  ungroup() %>% 
  select(id, year, time, end_dementia_death, everything())

data_long_dem %<>%
  group_by(id) %>%
  mutate(
    outcome_plr = ifelse(end_dementia_death == year & dementia == 1, 1, 0),
    competing_plr = ifelse(!is.na(dementia_date), 0, death_v)) %>%
  ungroup() %>%
  filter(!time == 0 & !competing_plr == 1) %>%
  arrange(id, time)

id2 <- data_long_dem %>% 
  filter(time == 0 & competing_plr == 1) %>% 
  pull(id)

data_long_dem %>% 
  count(outcome_plr, competing_plr)

### There are some differences between the wide structure and the long,
### There are 73 less cases of dementia because they had it after 21 yof
### but it is that a few had over 20 years of follow up

dem_a <- data_long_dem %>% filter(outcome_plr == 1) %>% pull(id)
dem_b <- data %>% filter(dementia == 1) %>% pull(id)

dif <- setdiff(dem_b, dem_a)

data %>% filter(id %in% dif) %>% View()
data %>% filter(id %in% dif) %>% count(t2dem_y >= 20)

data_long_dem %>% 
  filter(id %in% dif) %>% 
  count(id) %>% 
  filter(n < 19)

export(data_long_dem, here::here("01_data", "clean_data", "dementia_long.RData"))

