data_long <-
  rio::import(here::here("01_data", "clean_data", "dementia_long.RData"))

data_wide <-
  data_long %>%
  group_by(id) %>% 
  slice(max(row_number())) %>% 
  ungroup()

# merge with baseline covariates ------------------------------------------

baseline <-
  rio::import(here::here("01_data", "clean_data", "wide_noltfu.RData")) %>% 
  select(id, contains("1"))

data_wide <- data_wide %>% 
  left_join(baseline, by = c("id", "death_2015", "end_fup_2015", "e1"))

rio::export(data_wide, here::here("01_data", "clean_data", "wide_after_truncation.RData"))

# ## merge with cancer data for type of cancer
# 
# cancer_data <-
#   rio::import(here::here("01_data", "raw_data", "Oncodataset 19022020 R.RData"))

## Total n
total_n <- dim(data_wide)[1]

## Counts

total_cancer <- data_wide %>% 
  count(cancer_v) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(cancer_count = paste0(prop, "%", " (n=", n,")")) %>% 
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
  mutate(total = paste0(prop, "%", " (n=", n,")")) %>% 
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
  mutate(temp = paste0(m, " (IQR:", l, "-", u,")")) %>% 
  pull(temp)

dem_total <- data_wide %>% count(outcome_plr) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(temp = paste0(prop, "%", " (n=", n,")")) %>%
  filter(outcome_plr == 1) %>% 
  pull(temp)
  
dem_time <- data_long %>% 
  filter(dem_v == 1 & lag(dem_v) == 0) %>% 
  summarize(m = median(age_scale),
            l = quantile(age_scale, 0.25),
            u = quantile(age_scale, 0.75)) %>% 
  mutate(temp = paste0(m, " (IQR:", l, "-", u,")")) %>% 
  pull(temp)


## Baseline covs

mean_age <- data_wide %>% 
  summarize(temp = round(mean(age_0),  2)) %>% 
  #           b = round(sd(age_0),  2)) %>%
  # mutate(temp = paste0(a, " (SD:", b, ")")) %>% 
  pull(temp)


women_prop <- data_wide %>%
  count(sex) %>% 
  mutate(prop = round(100*n/sum(n),0)) %>% 
  mutate(total = paste0(prop, "%", " (n=", n,")")) %>% 
  filter(sex == 1) %>% 
  pull(total)

## Tableone

myvars <- c("sex", "age_0", "education", "apoe4", "smoke1", "bmi1", "oh1", "sbp1", 
          "ht1","hd_prev", "hd_v", "diabetes_prev", "diab_v",
          "stroke_prev", "stroke_v", "cancer_v") 

data_tableone <- data_wide %>% 
  select(myvars)

data_tableone %<>% 
  mutate(sex = ifelse(sex == 0, "Male", "Female"),
         education = case_when(
           education == 0 ~ "Lower",
           education == 1 ~ "Intermediate",
           education == 2 ~ "Higher",
           TRUE ~ "Unknown"),
         apoe4 = case_when(
           apoe4 == 0 ~ "Not carrier",
           apoe4 == 1 ~ "One allele carrier",
           apoe4 == 2 ~ "Two allele carrier"),
         smoke1 = case_when(
           smoke1 == 0 ~ "Never",
           smoke1 == 1 ~ "Former",
           smoke1 == 2 ~ "Current"),
         ht1 = ifelse(ht1 == 1, "History of hypertension", "No history of hypertension"),
         ht1 = as_factor(ht1),
         hd_prev = case_when(hd_prev == 1 ~ "History of heart disease",
                             hd_prev == 0 ~ "No history of heart disease"),
         hd_prev = as_factor(hd_prev),
         hd_v = case_when(hd_v == 1 ~ "Incident heart disease",
                          hd_v == 0 ~ "No incident heart disease"),
         hd_v = as_factor(hd_v),
         diabetes_prev = case_when(diabetes_prev == 1 ~ "History of diabetes",
                                   diabetes_prev == 0 ~ "No history of diabetes",
                                   diabetes_prev == 2 ~ "Unknown"),
         diabetes_prev = as_factor(diabetes_prev),
         diab_v = ifelse(diab_v == 1, "Incident diabetes", "No incident diabetes"),
         diab_v = fct_rev(as_factor(diab_v)),
         stroke_prev = ifelse(stroke_prev == 1, "History of stroke", "No history of stroke"),
         stroke_prev = as_factor(stroke_prev),
         stroke_v = ifelse(stroke_v == 1, "Incident stroke", "No incident stroke"),
         stroke_v = as_factor(stroke_v),
         cancer_v = ifelse(cancer_v == 1, "Incident cancer", "No incident cancer"),
         cancer_v = as_factor(cancer_v))



num <- c("age_0","bmi1", "oh1", "sbp1" )

cat <- myvars[!myvars %in% num]

library(tableone)
tableone <- CreateTableOne(vars = myvars, data = data_tableone, factorVars = cat)


tableone_exp <-
  print(
    tableone,
    explain = TRUE,
    test = FALSE,
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
  )

Characteristics <- rownames(tableone_exp)

tableone_complete <- tableone_exp %>% as_tibble() %>%
  cbind(Characteristics) %>% select(Characteristics, everything()) %>%
  mutate(
    Characteristics = str_replace(Characteristics, "smoke1", "Smoking status"),
    Characteristics = str_replace(Characteristics, "sex = ", ""),
    Characteristics = str_replace(Characteristics, "age_0", "Age at baseline"),
    Characteristics = str_replace(Characteristics, "education", "Educational attainment"),
    Characteristics = str_replace(Characteristics, "apoe4", "ApoE4"),
    Characteristics = str_replace(Characteristics, "bmi1", "Body Mass Index"),
    Characteristics = str_replace(Characteristics, "sbp1", "Systolic blood pressure (mmHg)"),
    Characteristics = str_replace(Characteristics, "ht1 = ", ""),
    Characteristics = str_replace(Characteristics, "hd_prev = ", ""),
    Characteristics = str_replace(Characteristics, "hd_v = ", ""),
    Characteristics = str_replace(Characteristics, "diab_v = ", ""),
    Characteristics = str_replace(Characteristics, "Unknown", "Unknown history of diabetes"),
    Characteristics = str_replace(Characteristics, "stroke_prev = ", ""),
    Characteristics = str_replace(Characteristics, "stroke_v = ", ""),
    Characteristics = str_replace(Characteristics, "cancer_v = ", ""),
  ) %>% 
  filter(!Characteristics %in% c("oh1 (mean (SD))", "diabetes_prev (%)")) %>% 
  filter(!Characteristics %in% str_detect(Characteristics, "No history of diabetes")) ## Not working :S

