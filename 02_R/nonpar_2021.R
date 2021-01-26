library(ggplot2)
library(tidyverse)
library(rio)
library(here)
library(wesanderson)

### Non parametric survival analysis

data_wide <- import(here("01_data", "no_ltfu", "wide_noltfu.RData"))

data_long_2014 <- import(here("01_data", "no_ltfu", "data_long_2014.RData"))

data_t2dem <- data_long_2014 %>% 
  mutate(year_dem_death = ifelse(!is.na(dementia_date),
                                 year(dementia_date),
                                 year(end_fup_2015))) %>% 
  filter(year <= year_dem_death) 

