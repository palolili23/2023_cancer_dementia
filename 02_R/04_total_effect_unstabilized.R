library(survival)
# library(survminer)
# library(broom)
library(dplyr)
library(splines)
library(magrittr)
library(here)
library(ggplot2)
# library(WeightIt)
# library(cobalt)

# Import data -------------------------------------------------------------
source(here::here("02_R", "04_auxiliary_fx.R"))

# data_long <-
#   load(here::here("01_data", "clean_data", "dementia_long.RData"))
# 
# data_wide <-
#   load(here::here("01_data", "clean_data", "wide_after_truncation.RData"))

data_long %<>% 
  group_by(id) %>% 
  mutate(
    fuptime = time,
    tstart = lag(fuptime),
    tstart = ifelse(is.na(tstart), -1, tstart)) %>% 
  ungroup()

data_long %<>% 
  mutate(outcome_plr = ifelse(competing_plr == 1, NA, outcome_plr))

data_wide %<>%
  mutate(dementia_efu =
           case_when(outcome_plr == 1 ~ 1,
                     competing_plr == 1 ~ 2,
                     TRUE ~ 0))

data_wide %<>%
  mutate(t2dem_efu =
           ifelse(t2dem >= 240, 240, t2dem))

data_wide %<>%
  mutate(combined_outcome_efu = ifelse(dementia_efu != 0, 1, 0))

data_long %<>% 
  mutate(both_outcomes = 
           case_when(outcome_plr == 1 ~ 1,
                     competing_plr == 1 ~2,
                     TRUE ~ 0))

data_long %<>% 
  mutate(combined_outcome = 
           case_when(outcome_plr == 1 ~ 1,
                     competing_plr == 1 ~ 1,
                     TRUE ~ 0))


# 1. Cancer Ever vs. Never ---------------------------------------------------


# 1.1. IPTW (cancer ever vs. never) ------------------------------------------

# By hand
cancer_den <-
  glm(
    cancer_v ~ bs(age_0) + sex + education + apoe4 + as.factor(smoke1),
    data = data_wide,
    family = binomial
  )

summary(cancer_den)

data_wide <- data_wide %>%
  mutate(
    p_denom = predict(cancer_den, type = "response"),
    w_cancer = ifelse(cancer_v == 1, 1/p_denom, 1/(1- p_denom)))

data_wide %<>%
  mutate(cancer_weight_t = ifelse(
    w_cancer > quantile(w_cancer, 0.99),
    quantile(w_cancer, 0.99),
    w_cancer
  ))

summary(data_wide$w_cancer)
summary(data_wide$cancer_weight_t)

# 1.2. IPCW for death, based on baseline covariates ----------------------------

## Death weights were once death = 1

death_den <-
  glm(
    competing_plr ~ bs(age_0, 3) + sex + education + apoe4 + 
      as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1,3) +
      as.factor(diabetes_prev) + cohort + cancer_v,
    data = data_wide,
    family = binomial
  )

data_wide %<>% mutate(
  p_death_den = predict(object = death_den, type = "response"),
  w_death = if_else(competing_plr == 0, 1/(1-p_death_den), 1))

data_wide %<>%
  mutate(w_death_t = ifelse(
    w_death > quantile(w_death, 0.99),
    quantile(w_death, 0.99),
    w_death
  ))


## Multiply weights

data_wide %<>%
  mutate(
    baseline_w = w_cancer* w_death,
    baseline_wt = cancer_weight_t * w_death_t)

# 1.a. AJ Unadjusted ------------------------------------------------------

cancer_ever_aj_crude <-survfit(
  Surv(
    t2dem_efu,
    event = as.factor(dementia_efu),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_wide)


rd_1a <- risks_km(cancer_ever_aj_crude) %>% 
  mutate(
    model = "Unadjusted"
  )

# 1.b. KM IPTW -----------------------------------------------

cancer_ever_aj <- survfit(
  Surv(
    t2dem_efu,
    event = as.factor(dementia_efu),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_wide,
  weights = cancer_weight_t,
  cluster = id)

rd_1b <- risks_km(cancer_ever_km) %>%
  mutate(model = "IPTW")


# 1.c. KM. IPCW for censoring -----------------------------------------------

cancer_ever_aj_ipcw <-survfit(
  Surv(
    t2dem_efu,
    event = as.factor(dementia_efu),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_wide,
  weights = baseline_wt,
  cluster = id)

rd_1c <- risks_km(cancer_ever_km_ipcw) %>%
  mutate(model = "IPTW + IPCW")


# 2. Time-varying cancer -----------------------------------------------------

# 2.1. Weights for t2cancer -----------------------------------------------

# We want the probability of having a cancer dx at each 
# time point for each participant who has not yet being diagnosed

data_long %<>% 
  group_by(id) %>%
  mutate(pastcancer = ifelse(cancer_v == 1 & lag(cancer_v) == 0, 0, cancer_v),
         pastcancer = ifelse(is.na(pastcancer), 0, pastcancer)) %>% 
  ungroup() 


# denominator
mod <- glm(cancer_v ~ bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
             as.factor(smoke_lag) + bs(sbp_lag, 3) + bs(bmi_lag, 3) + 
             ht_lag + ht_drug_lag + diab_lag + cohort,
           family = quasibinomial(), data = subset(data_long , 
                                                   pastcancer == 0))

summary(mod)

data_long %<>%
  mutate(pred_a_den = ifelse(pastcancer == 1, 
                             1, # the pr of ever having cancer is 1 for those with cancer
                             predict(mod, type = 'response'))) 

data_long %<>% 
  mutate(den_a = ifelse(cancer_v == 1, pred_a_den, 1 - pred_a_den)) %>% 
  group_by(id) %>% 
  mutate(dencum = cumprod(den_a)) %>% 
  ungroup() %>% 
  mutate(w_cancer = 1/dencum)


data_long %<>%
  mutate(w_cancer_t = ifelse((w_cancer > quantile(w_cancer, 0.99)),
                             quantile(w_cancer, 0.99),
                             w_cancer))
# 
# data_long %>% 
#   ggplot(aes(x = as_factor(time), y = sw_cancer_t)) +
#   geom_boxplot()

# 2.2. IPCW. Weights on death - dementia ---------------------------------------------

death_den <- glm(
  competing_plr ~ cancer_lag + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
    as.factor(smoke_lag) + bs(sbp_lag, 3) + bs(bmi_lag, 3) + 
    ht_lag + ht_drug_lag + hd_lag + stroke_lag + diab_lag + cohort,
  data = data_long,
  family = quasibinomial
)

broom::tidy(death_den, exponentiate = TRUE)

data_long$p_denom = predict(death_den, data_long, type = "response")

data_long %<>%
  group_by(id) %>%
  mutate(
    w_death = 1 / cumprod(1 - p_denom)    # 1 - p because is the probability of not dying
  ) %>%
  ungroup()

data_long %<>%
  mutate(
    w_death_t = ifelse((w_death > quantile(w_death, 0.99)),
                       quantile(w_death, 0.99), w_death))

data_long %>% 
  filter(!is.na(outcome_plr)) %>%  
  ggplot(aes(x = as_factor(time), y = sw_death_t)) +
  geom_boxplot()

data_long %<>%
  mutate(weights_both = w_cancer * w_death, 
         weights_both_t = w_cancer_t * w_death_t)

data_long %>% filter(tstart >= fuptime) %>% select(id, time, tstart, fuptime)

# 3.a. KM unadjusted ------------------------------------------------------

km_t2cancer <- survfit(
  Surv(
    tstart,
    time2 = fuptime,
    event = as.factor(both_outcomes),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_long,
  cluster = id,
  id = id
)

rd_3a <- risks_km(km_t2cancer) %>% 
  mutate(model = "Unadjusted")


# 3.b. KM IPTW 2tcancer  ------------------------------

aj_t2cancer_iptw <- 
  km_t2cancer <- survfit(
    Surv(
      tstart,
      time2 = fuptime,
      event = as.factor(both_outcomes),
      type = 'mstate'
    ) ~ cancer_v,
    data = data_long,
    cluster = id,
    id = id, 
weights = w_cancer_t)

plot_3b <- plot_km(aj_t2cancer_iptw, "Risk of dementia, time to cancer") + 
  labs(subtitle = "IPTW")

rd_3b <- risks_km(km_t2cancer_iptw) %>% 
  mutate(model = "IPTW")

# 3.c. KM IPCW T2cancer  ---------------

km_t2cancer_ipcw <- 
  aj_t2cancer_iptw <- 
  km_t2cancer <- survfit(
    Surv(
      tstart,
      time2 = fuptime,
      event = as.factor(both_outcomes),
      type = 'mstate'
    ) ~ cancer_v,
    data = data_long,
    cluster = id,
    id = id, weights = weights_both_t)



results_output_unstabilized <- list(cancer_ever_km_crude, 
                                    cancer_ever_km, 
                                    km_t2cancer, km_t2cancer_iptw, 
                                    km_t2cancer_ipcw)

