## Plots 2021
library(ggplot2)
library(tidyverse)
library(rio)
library(wesanderson)
library(magrittr)
palette <- wes_palette("Moonrise3")

## Import dataset that keeps deads until 2014

data_wide <- import(here::here("01_data", "clean_data", "wide_noltfu.RData"))

data_long <- import(here::here("01_data", "clean_data", "data_long_2015.RData"))
data_dem <- import(here::here("01_data", "clean_data", "dementia_long.RData"))

# data_long %>% select(age_0, age_scale) %>% View()

# data_long %>%
#   filter(id == 17002) %>%
#   select(
#     id,
#     cohort,
#     year,
#     e1,
#     cancer,
#     cancer_v,
#     cancer_date,
#     dementia,
#     dem_v,
#     dem_date,
#     death_2015,
#     death_v,
#     death_year,
#     status,
#     t2death_y
#   ) %>% View()

## Plot distribution of events over time

## Keeps those who died
data_long %>% 
  group_by(age_scale) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n)),0) %>% 
  ungroup() %>% 
  filter(age_scale <= 90) %>% 
  ggplot(aes(age_scale, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  theme_minimal()

data_wide %>% 
  ggplot(aes(age_0)) +
  geom_histogram()

## Removes those who died after year of death
data_long %>% 
  filter(year <= year(end_fup_2015)) %>%
  group_by(age_scale) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n))) %>% 
  ungroup() %>%
  filter(age_scale <= 90) %>% 
  ggplot(aes(age_scale, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  theme_minimal()

library(lubridate)

## Removes those who had dementia after year of dementia dx
data_dem %>% 
  group_by(age_scale) %>% 
  count(status) %>%
  mutate(prop = round(100*n/sum(n))) %>% 
  ungroup() %>% 
  ggplot(aes(age_scale, prop, fill = status)) +
  geom_area() + 
  scale_fill_manual(values = palette) +
  theme_minimal()

##### Plot longitudinal distribution over time

data_id <- data_long %>%ungroup() %>% count(id) %>% sample_n(100) %>% pull(id)

data_long %>%
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


