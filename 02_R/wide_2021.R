library(rio)
library(lubridate)
library(tidyverse)
library(here)
library(magrittr)
# library(summarytools)

# Import first cohort -----------------------------------------------------

data <- import(here::here("01_data", "clean_data.RData"))

data <- data %>%  
  filter(!is.na(e1)) %>% 
  filter(dementia_inc != 2) %>% 
  filter(dementia_at_risk == 1) %>% 
  #filter(cancer_prev != 1) %>% # KIM this one is not needed anymore
  filter(!is.na(education_three_levels)) %>% 
  filter(!is.na(smoke1)) %>% 
  filter(!is.na(bmi1)) %>% 
  filter(!is.na(sbp1)) %>% 
  filter(!is.na(ht_drug1)) %>% 
  filter(!is.na(censor_date)) %>% 
  rename(id = ergoid, cohort = rs_cohort) 

#KIM: SELECT ONLY THOSE IN RS-1 AND RS-2
data <- subset(data,data$cohort!=3)

#KIM: SELECT ONLY THOSE WITH IC BASED ON STROKE DATASET
data <- subset(data,data$ic_ok_2016==1)


## truncate death in 2015
data <- data %>%
  mutate(
    death = ifelse(!is.na(mort_date), 1, 0),
    end_fup_2015 = as_date(ifelse(
      censor_date > ymd("2014-12-31"),
      ymd("2014-12-31"),
      censor_date
    )),
    death_2015 = ifelse(mort_date > ymd("2014-12-31") &
                          !is.na(mort_date), 0, death)
  )

### Calculate t2events in days, months and years

data %<>%
  mutate(
    t2death = round(as.numeric(as.period((e1 %--% end_fup_2015), "months"
    ), "months"), 0),
    t2dem = round(as.numeric(as.period((e1 %--% dementia_date), "months"
    ), "months"), 0),
    t2cancer = round(as.numeric(as.period((e1 %--% cancer_date), "months"
    ), "months"), 0))

data %<>%
  mutate(
    t2death_y = round(as.numeric(as.period((e1 %--% end_fup_2015), "years"
    ), "years"), 3),
    t2dem_y = round(as.numeric(as.period((e1 %--% dementia_date), "years"
    ), "years"), 3),
    t2cancer_y = round(as.numeric(as.period((e1 %--% cancer_date), "years"
    ), "years"), 3))

data %<>%
  mutate(
    t2death_d = round(as.numeric(as.period((e1 %--% end_fup_2015), "days"
    ), "days"), 0),
    t2dem_d = round(as.numeric(as.period((e1 %--% dementia_date), "days"
    ), "days"), 0),
    t2cancer_d = round(as.numeric(as.period((e1 %--% cancer_date), "days"
    ), "days"), 0))


## t2dem and t2cancer are empty if date variables are empty, so use the end_fup_2015 to fill the fup time

## dementia
data <- data %>%
  mutate(
    dementia = case_when(
      !is.na(dementia_date) & dementia_date <= end_fup_2015 ~ 1,
      death_2015 == 1  &
        (is.na(dementia_date) | dementia_date > end_fup_2015) ~ 2,
      TRUE ~ 0
    ),
    t2dem = ifelse((is.na(dementia_date) | dementia_date > end_fup_2015), t2death, t2dem),
    t2dem_y = ifelse((is.na(dementia_date) |dementia_date > end_fup_2015), t2death_y, t2dem_y),
    t2dem_d = ifelse((is.na(dementia_date) |dementia_date > end_fup_2015), t2death_d, t2dem_d),
  ) 


## Cancer
data <- data %>%
  mutate(
    cancer = case_when(!is.na(cancer_date) & cancer_date <= end_fup_2015 ~ 1,
      TRUE ~ 0),
    t2cancer = ifelse((is.na(cancer_date) | cancer_date > end_fup_2015), t2death, t2cancer),
    t2cancer_y = ifelse((is.na(cancer_date) |cancer_date > end_fup_2015), t2death_y, t2cancer_y),
    t2cancer_d = ifelse((is.na(cancer_date) | cancer_date > end_fup_2015), t2death_d, t2cancer_d))


data %>% 
  select(cancer, dementia, death, death_2015, t2dem, t2cancer, t2death) %>% View()

### There are some individuals with negative time, all of them have prevalent cancer
### Two individuals had cancer_prev = 1, cancer_inc = 1, and t2cancer and t2 death have the same date
### Until this point, these and further (dem dx at t1, or dem dx and death dx at same time) will not be addressed


# data$t2cancer <- ifelse(data$t2cancer<0&data$cancer_prev==1,0,data$t2cancer)
# 
# data$t2cancer_y <- ifelse(data$t2cancer_y<0&data$cancer_prev==1,0,data$t2cancer_y)


prev_cancer <- data %>% 
  filter(cancer_prev == 1 | t2cancer < 0) %>% 
  select(id, cancer_prev, cancer_inc, cancer, death_2015, dementia, 
         t2cancer, t2dem, t2death)

neg_t2cancer <- data %>%
  filter(t2cancer < 0) %>% 
  select(id, cancer_prev, cancer_inc,cancer, death_2015, dementia, 
         t2cancer, t2dem, t2death)

prev_cancer %>% setdiff(neg_t2cancer)

export(data, here::here("01_data", "no_ltfu", "wide_noltfu.Rdata"))

data %>% 
  filter(t2dem_y <1) %>% select(cancer_prev, cancer_inc,cancer, death_2015, dementia, 
                                t2cancer, t2dem, t2death)
