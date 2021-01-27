## Plots 2021
library(ggplot2)
library(tidyverse)
library(rio)
library(wesanderson)
library(magrittr)
palette <- wes_palette("Moonrise3")

## Import dataset that keeps deads until 2014

data_wide <- import(here::here("01_data", "clean_data", "wide_noltfu.RData"))

data_long_2014 <- import(here::here("01_data", "clean_data", "data_long_2014.RData"))

data_long_2014 %<>%
  mutate(age_scale = age_0 + (row_number() - 1))

data_long_2014 %>% 
  filter(id == 17002) %>% 
  select(id, cohort, year, e1, cancer, cancer_v, cancer_date, dementia, dem_v, dem_date, death_2015, death_v, death_year, status, t2death_y)%>% View() 

## Plot distribution of events over time

## Keeps those who died
data_long_2014 %>% 
  # filter(between(age_0, 60, 75)) %>%
  group_by(cohort, year) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n))) %>% 
  ungroup() %>% 
  filter(cohort == 1) %>% 
  ggplot(aes(year, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  # facet_wrap( . ~ cohort, scales = "free") +
  theme_minimal()

data_wide %>% 
  ggplot(aes(age_0)) +
  geom_histogram() +
  facet_wrap(.~cohort)

## Removes those who died after year of death
data_long_2014 %>% 
  filter(between(age_0, 60, 75)) %>%
  # filter(year <= year(end_fup_2015)) %>% 
  group_by(cohort, year) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n))) %>% 
  ungroup() %>% 
  ggplot(aes(year, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  facet_wrap(.~ cohort, scales = "free") +
  theme_minimal()

library(lubridate)
## Removes those who had dementia after year of dementia dx
data_long_2014 %>% 
  mutate(year_dem_death = ifelse(!is.na(dementia_date),
                                 year(dementia_date),
                                 year(end_fup_2015))) %>% 
  filter(year <= year_dem_death) %>% 
  filter(year <= 2014) %>% 
  # filter(between(age_0, 60, 75)) %>%
  # filter(!status %in% c("Alive, free of cancer and dementia", "Dead")) %>% 
  group_by(cohort, year) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n))) %>% 
  ungroup() %>% 
  ggplot(aes(year, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  facet_wrap(.~ cohort, scales = "free") +
  theme_minimal()

##### Plot longitudinal distribution over time

data_id <- data_long_2014 %>%ungroup() %>% count(id) %>% sample_n(100) %>% pull(id)

data_long_2014 %>%
  group_by(id) %>% 
  mutate(age_death = max(age_scale)) %>% 
  ungroup() %>% 
  # filter(id %in% data_id) %>% 
  filter(cohort == 1) %>% 
  ggplot(aes(
    x = age_scale,
    y = fct_reorder(as.factor(id), age_death),
    group = id,
   color = as.factor(status))) +
  geom_line() +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_wrap(.~cancer)


