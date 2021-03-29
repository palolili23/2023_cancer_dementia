library(rio)
library(lubridate)
library(tidyverse)
library(here)
library(magrittr)
# library(summarytools)

# Import first cohort -----------------------------------------------------

data <- import(here::here("01_data", "clean_data","complete_data.RData"))

data <- data %>%  
  rename(id = ergoid, cohort = rs_cohort) %>% 
  filter(!is.na(e1)) %>% 
  filter(dementia_inc != 2) %>% 
  filter(dementia_at_risk == 1) %>% 
  filter(cancer_prev != 1) %>% # Based on the ideal trial
  # filter(!is.na(education_three_levels)) %>% 
  # filter(!is.na(smoke1)) %>% 
  filter(!is.na(bmi1)) %>%
  filter(!is.na(sbp1)) %>%
  filter(!is.na(ht_drug1)) %>%
  filter(!is.na(censor_date)) %>% 
  filter(between(age_0, 60, 70)) %>% 
  filter(cohort == 1)

#KIM: SELECT ONLY THOSE WITH IC BASED ON STROKE DATASET
data <- subset(data,data$ic_ok_2016==1)

data %>% ggplot(aes(age_0)) + geom_histogram()

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

data %>% 
  group_by(cohort) %>% 
  summarize(
    max(t2death_y),
    max(t2dem_y, na.rm = TRUE),
    max(t2cancer_y, na.rm = TRUE))

## t2dem and t2cancer are empty if date variables are empty,
## so use the end_fup_2015 to fill the fup time

## dementia
data <- data %>%
  mutate(
    dementia = case_when(
      !is.na(dementia_date) & dementia_date <= end_fup_2015 ~ 1,
      death_2015 == 1  &
        (is.na(dementia_date) | dementia_date > end_fup_2015) ~ 2,
      TRUE ~ 0
    ),
    t2dem = ifelse((is.na(dementia_date) | dementia_date > end_fup_2015), 
                   t2death, t2dem),
    t2dem_y = ifelse((is.na(dementia_date) |dementia_date > end_fup_2015), 
                     t2death_y, t2dem_y),
    t2dem_d = ifelse((is.na(dementia_date) |dementia_date > end_fup_2015), 
                     t2death_d, t2dem_d),
  ) 


## Cancer
data <- data %>%
  mutate(
    cancer = ifelse(!is.na(cancer_date) & cancer_date <= end_fup_2015, 1, 0),
    t2cancer = ifelse((is.na(cancer_date) | cancer_date > end_fup_2015), 
                      t2death, t2cancer),
    t2cancer_y = ifelse((is.na(cancer_date) |cancer_date > end_fup_2015), 
                        t2death_y, t2cancer_y),
    t2cancer_d = ifelse((is.na(cancer_date) | cancer_date > end_fup_2015), 
                        t2death_d, t2cancer_d))


# data %>% 
#   select(cancer, dementia, death, death_2015, t2dem, t2cancer, t2death,
#          cancer_date, dementia_date, end_fup_2015, mort_date, e1) %>% View()

### Code check-----------------------------------------------------------

### There are some individuals with negative time,
### all of them have prevalent cancer
### Two individuals had cancer_prev = 1, cancer_inc = 1,
### and t2cancer and t2 death have the same date
### Until this point, these and further (dem dx at t1, or dem dx
### and death dx at same time) will not be addressed

# prev_cancer <- data %>%
#   filter(cancer_prev == 1 | t2cancer < 0) %>%
#   select(id, cancer_prev, cancer_inc, cancer, death_2015, dementia,
#          t2cancer, t2dem, t2death)
# 
# neg_t2cancer <- data %>%
#   filter(t2cancer < 0) %>%
#   select(id, cancer_prev, cancer_inc,cancer, death_2015, dementia,
#          t2cancer, t2dem, t2death)
# 
# prev_cancer %>% setdiff(neg_t2cancer)


# Make categories for missing data ----------------------------------------

### Missing data

data %>% count(is.na(education_three_levels))
data %>% count(is.na(apoe4))
data %>% count(is.na(ht_drug1))
data %>% count(is.na(hd_prev))
data %>% count(cancer_prev)
data %>% count(ht1)

data %<>% 
  rename(education = education_three_levels) %>% 
  mutate(
    education = ifelse(is.na(education), "Unknown", education),
    apoe4 = ifelse(is.na(apoe4),"Unknown" , apoe4),
    ht_drug1 = ifelse(is.na(ht_drug1), 2, ht_drug1),
    hd_prev = ifelse(is.na(hd_prev), 2, hd_prev),)

# Wide truncated 20 years -------------------------------------------------

data %>% summarise(max(t2dem))
data %>% summarise(max(t2cancer))
data %>% summarise(max(t2death))

data %>% count(cancer, dementia)
# 
# data_test <- data %>% 
#   mutate(
#     dementia_20 = ifelse(t2dem > 240, 0, dementia),
#     death_20 = ifelse(t2death > 240, 0, death_2015),
#     cancer_20 = ifelse(t2cancer > 240, 0, cancer))
# 
# data %>% count(cancer, dementia_20)
# data %>% count(cancer_20, dementia_20)
# 
# data <- data %>% 
#   mutate(
#     t2dem_20 = ifelse(t2dem > 240, 240, t2dem),
#     t2death_20 = ifelse(t2death > 240, 240, t2death),
#     t2cancer_20 = ifelse(t2cancer > 240, 240, t2death))

### Those who have cancer after dementia need to be changed

# data %>% filter(t2dem < t2cancer & !is.na(cancer_date) & 
#                   !is.na(dementia_date)) %>%
#   select(id,
#          starts_with("t2"),
#          contains("cancer"),
#          contains("death"),
#          contains("dem")) %>%
#   view()

data <- data %>%
  mutate(
    cancer = ifelse(
      (t2dem < t2cancer) & !is.na(cancer_date) & !is.na(dementia_date),
      0,
      cancer
    ),
    t2cancer = ifelse(
      (t2dem < t2cancer) &
        !is.na(cancer_date) & !is.na(dementia_date),
      t2dem,
      t2cancer
    )
  )


export(data, here::here("01_data", "clean_data","wide_noltfu.Rdata"))

# Export ------------------------------------------------------------------


data %>% 
  filter(t2dem <1) %>% 
  select(cancer_prev, cancer_inc,cancer, death_2015,
                              dementia, t2cancer, t2dem, t2death)
