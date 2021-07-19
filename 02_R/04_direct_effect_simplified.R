library(tidyverse)
library(rio)
library(survival)
library(survminer)
library(broom)
library(splines)
library(magrittr)
library(WeightIt)
library(cobalt)

# Import data -------------------------------------------------------------
source(here::here("02_R", "04_auxiliary_fx.R"))

data_long <-
  import(here::here("01_data", "clean_data", "dementia_long.RData"))

data_wide <-
  import(here::here("01_data", "clean_data", "wide_after_truncation.RData"))

data_long %<>% 
  group_by(id) %>% 
  mutate(
    fuptime = time,
    tstart = lag(fuptime),
    tstart = ifelse(is.na(tstart), -1, tstart)) %>% 
  ungroup()

## idea of how numbers look like over follow-up

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

## If competing_plr = 1, then outcome_plr = NA

data_long %<>% 
  mutate(outcome_plr = ifelse(competing_plr == 1, NA, outcome_plr))

# 1. Cancer Ever vs. Never ---------------------------------------------------


# 1.1. IPTW (cancer ever vs. never) ------------------------------------------

# By hand
# cancer_den <-
#   glm(
#     cancer_v ~ bs(age_0) + sex + education + apoe4 + as.factor(smoke1),
#     data = data_wide,
#     family = binomial
#   )
# 
# summary(cancer_den)
# 
# cancer_num <- glm(cancer_v ~ 1, data = data_wide)
# 
# summary(cancer_num)
# 
# data_wide <- data_wide %>%
#   mutate(
#     p_num = predict(cancer_num, type = "response"),
#     p_denom = predict(cancer_den, type = "response"),
#     w_cancer = ifelse(cancer_v == 1, p_num/p_denom, (1 - p_num)/(1- p_denom)))


## Weightit package, makes nice loveplots

w.out1 <-
  weightit(
    cancer_v ~ bs(age_0, 3) + sex + education + as.factor(smoke1) + cohort,
    data = data_wide,
    stabilize = TRUE,
    estimand = "ATE",
    method = "ps"
  )

w.out1

data_wide <- data_wide %>% 
  mutate(cancer_weight= w.out1$weights)

love.plot(w.out1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#011A5E", "#e4a803"))

data_wide %<>%
  mutate(cancer_weight_t = ifelse(
    cancer_weight > quantile(cancer_weight, 0.99),
    quantile(cancer_weight, 0.99),
    cancer_weight
  ))

summary(data_wide$cancer_weight)
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

death_num <-
  glm(
    competing_plr ~ 1,
    data = data_wide,
    family = binomial
  )

data_wide %<>% mutate(
  p_death_den = predict(object = death_den, type = "response"),
  p_death_num = predict(object = death_num, type = "response"),
  sw_death = if_else(competing_plr == 0, (1-p_death_num)/(1-p_death_den), 1))

data_wide %<>%
  mutate(sw_death_t = ifelse(
    sw_death > quantile(sw_death, 0.99),
    quantile(sw_death, 0.99),
    sw_death
  ))


## Multiply weights

data_wide %<>%
  mutate(
    baseline_w = cancer_weight * sw_death,
    baseline_wt = cancer_weight_t * sw_death_t)

# 1.a. KM Unadjusted ------------------------------------------------------

cancer_ever_km_crude <-
  survfit(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide)

plot_km(cancer_ever_km_crude, "KM crude") +
  labs(subtitle = "Under unconditional independent censoring of death")

cancer_ever_cox_crude <-
  coxph(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide)

hr_1a <- cancer_ever_cox_crude %>% 
  tidy_hr() %>% 
  mutate(
    model = "Unadjusted"
  )

rd_1a <- risks_km(cancer_ever_km_crude) %>% 
  mutate(
    model = "Unadjusted"
  )


# 1.b. KM IPTW -----------------------------------------------

cancer_ever_km <-
  survfit(
    Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
    data_wide,
    weights = cancer_weight_t,
    cluster = id
  )

plot_1b <- plot_km(cancer_ever_km, "Risk of dementia, ever vs. never") +
  labs(subtitle = "IPTW")

cancer_ever_cox <-
  coxph(
    Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
    data_wide,
    weights = cancer_weight_t,
    cluster = id
  )

hr_1b <- cancer_ever_cox %>%
  tidy_hr() %>% mutate(model = "IPTW")

rd_1b <- risks_km(cancer_ever_km) %>%
  mutate(model = "IPTW")

# 1.c. KM. IPCW for censoring -----------------------------------------------

cancer_ever_km_ipcw <-
  survfit(
    Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
    data_wide,
    weights = baseline_wt,
    cluster = id
  )

plot_1c <- plot_km(cancer_ever_km_ipcw, "Risk of dementia, ever vs. never") +
  labs(subtitle = "IPTW + IPCW")

cancer_ever_cox_ipcw <-
  coxph(
    Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
    data_wide,
    weights = baseline_wt,
    cluster = id
  )

hr_1c <- cancer_ever_cox_ipcw %>%
  tidy_hr() %>%
  mutate(model = "IPTW + IPCW")

rd_1c <- risks_km(cancer_ever_km_ipcw) %>%
  mutate(model = "IPTW + IPCW")

hr_ever_never <- bind_rows(hr_1a, hr_1b, hr_1c) %>% 
  mutate(Proxy = "Ever vs. Never") 


# 1.d. Lower bound --------------------------------------------------------

km_ever_ipcw_low_bound <- survfit(
  Surv(
    t2dem_efu,
    event = as.factor(dementia_efu),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_wide,
  weights = cancer_weight_t)

rd_1d <- lower_bound(km_ever_ipcw_low_bound)

# 1.e. Upper bound -------------------------------------------------------------

km_ever_ipcw_up_bound <- survfit(
  Surv(
    t2dem_efu,
    event = combined_outcome_efu,
  ) ~ cancer_v,
  data = data_wide,
  weights = cancer_weight_t)

rd_1e  <- km_ever_ipcw_up_bound %>% 
  risks_km()  %>% 
  mutate(model = "Upper bound")


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

# numerator

mod_num <- glm(cancer_v ~ bs(time, 3),
               family = quasibinomial(),
               data = subset(data_long, pastcancer == 0))

data_long %<>% 
  mutate(pred_a_num = ifelse(pastcancer == 1, 
                             1, # the pr of ever having cancer is 1 for those with cancer
                             predict(mod_num, type = 'response'))) 

data_long %<>% 
  mutate(num_a = ifelse(cancer_v == 1, pred_a_num, 1 - pred_a_num),
         den_a = ifelse(cancer_v == 1, pred_a_den, 1 - pred_a_den)) %>% 
  group_by(id) %>% 
  mutate(numcum = cumprod(num_a),
         dencum = cumprod(den_a)) %>% 
  ungroup() %>% 
  mutate(sw_cancer = numcum/dencum,
         w_cancer = 1/dencum)


data_long %<>%
  mutate(
    sw_cancer_t = ifelse((sw_cancer > quantile(sw_cancer, 0.99)),
                           quantile(sw_cancer, 0.99),
                         sw_cancer
    ),
    w_cancer_t = ifelse((w_cancer > quantile(w_cancer, 0.99)),
                          quantile(w_cancer, 0.99),
                        w_cancer
    )
  )

data_long %>% 
  ggplot(aes(x = as_factor(time), y = sw_cancer_t)) +
  geom_boxplot()

# 2.2. IPCW. Weights on death - dementia ---------------------------------------------

death_den <- glm(
  competing_plr ~ cancer_lag + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
    as.factor(smoke_lag) + bs(sbp_lag, 3) + bs(bmi_lag, 3) + 
    ht_lag + ht_drug_lag + hd_lag + stroke_lag + diab_lag + cohort,
  data = data_long,
  family = quasibinomial
)

broom::tidy(death_den, exponentiate = TRUE)

death_num <-
  glm(
    competing_plr ~ cancer_lag + bs(time, 3) + cohort,
    data = data_long,
    family = quasibinomial
  )

broom::tidy(death_num, exponentiate = TRUE)

data_long$p_denom = predict(death_den, data_long, type = "response")

data_long$p_num = predict(death_num, data_long, type = "response")

data_long %<>%
  group_by(id) %>%
  mutate(
    w_death = 1 / cumprod(1 - p_denom),    # 1 - p because is the probability of not dying
    sw_death = cumprod(1 - p_num) / cumprod(1 - p_denom)
  ) %>%
  ungroup()

data_long %<>%
  mutate(
    w_death_t = ifelse((w_death > quantile(w_death, 0.99)),
                       quantile(w_death, 0.99), w_death),
    sw_death_t = ifelse((sw_death > quantile(sw_death, 0.99)),
                        quantile(sw_death, 0.99), sw_death)
  )

data_long %>% 
  filter(!is.na(outcome_plr)) %>%  
  ggplot(aes(x = as_factor(time), y = sw_death_t)) +
  geom_boxplot()

data_long %<>%
  mutate(weights_both = w_cancer * w_death, 
         weights_both_t = w_cancer_t * w_death_t,
         sw_weights_both = sw_cancer * sw_death,
         sw_weights_both_t = sw_cancer_t * sw_death_t)

data_long %>% filter(tstart >= fuptime) %>% select(id, time, tstart, fuptime)

# 3.a. KM unadjusted ------------------------------------------------------

km_t2cancer <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id)

plot_3a <- plot_km(km_t2cancer, "Risk of dementia, time to cancer") + 
  labs(subtitle = "Unadjusted")


cox_t2cancer <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id)

hr_3a <- cox_t2cancer %>% 
  tidy_hr() %>% 
  mutate(model = "Unadjusted")

rd_3a <- risks_km(km_t2cancer) %>% 
  mutate(model = "Unadjusted")

# 3.b. KM IPTW 2tcancer  ------------------------------

km_t2cancer_iptw <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_cancer_t)

plot_3b <- plot_km(km_t2cancer_iptw, "Risk of dementia, time to cancer") + 
  labs(subtitle = "IPTW")


cox_t2cancer_iptw <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_cancer_t)

hr_3b <- cox_t2cancer_iptw %>% 
  tidy_hr() %>% 
  mutate(model = "IPTW")

rd_3b <- risks_km(km_t2cancer_iptw) %>% 
  mutate(model = "IPTW")

# 3.c. KM IPCW T2cancer  ---------------

km_t2cancer_ipcw <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_weights_both_t)

cox_t2cancer_ipcw <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_weights_both_t)

plot_3c <- plot_km(km_t2cancer_ipcw, "Risk of dementia, time to cancer") + 
  labs(subtitle = "IPTW + IPCW")

hr_3c <- cox_t2cancer_ipcw %>% 
  tidy_hr() %>% 
  mutate(model = "IPTW + IPCW")

rd_3c <- risks_km(km_t2cancer_ipcw) %>% 
  mutate(model = "IPTW + IPCW")

hr_t2cancer <- bind_rows(hr_3a, hr_3b, hr_3c)%>% 
  mutate(Proxy = "Time to cancer") 


# 3.d. Lower bound ------------------------------------------------------------------

km_t2c_ipcw_low_bound <- survfit(
  Surv(
    tstart,
    time2 = fuptime,
    event = as.factor(both_outcomes),
    type = 'mstate'
  ) ~ cancer_v,
  data = data_long,
  cluster = id,
  id = id,
  weights = sw_cancer_t
)

rd_3d  <- lower_bound(km_t2c_ipcw_low_bound)

# 3.e. Upper bound -------------------------------------------------------------

km_t2cancer_ipcw_up_bound <- survfit(
  Surv(
    tstart,
    time2 = fuptime,
    event = combined_outcome,
  ) ~ cancer_v,
  data = data_long,
  cluster = id,
  id = id,
  weights = sw_cancer_t
)

rd_3e  <- km_t2cancer_ipcw_up_bound %>% 
  risks_km()  %>% 
  mutate(model = "Upper bound")



# 4. All results ----------------------------------------------------------

hr_results <- bind_rows(hr_ever_never, hr_t2cancer)

export(hr_results, here::here("02_R", "hr_results.csv"))


## Bounds

bounds_ever_cancer <- bind_rows(rd_1c, rd_1d, rd_1e) %>% 
  select(model, everything()) %>% 
  arrange(rd) %>% 
  mutate(Proxy = "Ever vs. never")

bounds_t2vcancer <- bind_rows(rd_3c, rd_3d, rd_3e) %>% 
  select(model, everything()) %>% 
  arrange(rd) %>% 
  mutate(Proxy = "Time to cancer")

table_results_bounds <- bind_rows(bounds_ever_cancer, bounds_t2vcancer) %>% 
  select(Proxy, everything())

export(table_results_bounds, here::here("02_R", "table_results_bounds.csv"))

plots <- list(
  plot_1b, plot_1c, plot_3b, plot_3c)

save(plots, file = here::here("03_figs", "plots.RData"))

models <- list(cancer_ever_km_ipcw, km_t2cancer_ipcw)

save(models, file = here::here("02_R", "models.RData"))
