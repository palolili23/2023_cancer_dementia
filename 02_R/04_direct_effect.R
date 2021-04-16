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
source(here::here("02_R", "04b_auxiliary_fx.R"))

data_long <-
  import(here::here("01_data", "clean_data", "dementia_long.RData"))

data_wide <-
  import(here::here("01_data", "clean_data", "wide_after_truncation.RData"))

data_long %<>% 
  group_by(id) %>% 
  mutate(
    fuptime = time - 1,
    tstart = lag(fuptime),
    tstart = ifelse(is.na(tstart), -1, tstart)) %>% 
  ungroup()

## idea of how numbers look like over follow-up

data_long %<>% 
  mutate(both_outcomes = 
           case_when(outcome_plr == 1 ~ 1,
                     competing_plr == 1 ~2,
                     TRUE ~ 0))

data_wide %<>%
  mutate(dementia_efu =
           case_when(outcome_plr == 1 ~ 1,
                     competing_plr == 1 ~ 2,
                     TRUE ~ 0))

data_wide %<>%
  mutate(t2dem_efu =
           ifelse(t2dem >= 240, 240, t2dem))


# 1. Cancer Ever vs. Never ---------------------------------------------------


# 1.1. IPTW (cancer ever vs. never) ------------------------------------------

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


## Check standardized mean

w.out1 <-
  weightit(
    cancer_v ~ bs(age_0, 3) + sex + education + apoe4 + as.factor(smoke1) +
      bs(sbp1, 3) + bs(bmi1, 3) + ht1 + as.factor(diabetes_prev) + cohort,
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

w.out_death <-
  weightit(
    competing_plr ~ bs(age_0, 3) + sex + education + apoe4 + 
      as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1,3) +
      as.factor(diabetes_prev) +  cancer_v + cohort,
    data = data_wide,
    stabilize = TRUE,
    estimand = "ATE",
    method = "ps"
  )

w.out_death

data_wide <- data_wide %>% 
  mutate(death_weight = w.out_death$weights)

love.plot(w.out_death) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#011A5E", "#e4a803"))

data_wide %<>%
  mutate(death_weight_t = ifelse(
    death_weight > quantile(death_weight, 0.99),
    quantile(death_weight, 0.99),
    death_weight
  ))


summary(data_wide$death_weight)
summary(data_wide$death_weight_t)

## Multiply weights

data_wide %<>%
  mutate(
    baseline_w = cancer_weight * death_weight,
    baseline_wt = cancer_weight_t * death_weight_t)

# 1.a. KM. Unadjusted for confounding

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


# 1.b. KM Independent censoring -----------------------------------------------

cancer_ever_km <-
  survfit(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide, weights = cancer_weight_t, cluster = id)

plot_km(cancer_ever_km, "Direct effect") +
  labs(subtitle = "Under unconditional independent censoring of death")

cancer_ever_cox <-
  coxph(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide, weights = cancer_weight_t, cluster = id)

hr_1b <- cancer_ever_cox %>% 
  tidy_hr() %>% mutate(
    model = "IPTW"
  )

rd_1b <- risks_km(cancer_ever_km) %>% 
  mutate(
    model = "IPTW"
  )

# 1.c. KM. IPCW for censoring -----------------------------------------------

cancer_ever_km_ipcw <-
  survfit(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide, weights = baseline_wt, cluster = id)

plot_km(cancer_ever_km_ipcw, "Direct effect") +
  labs(subtitle = "With IPCW for death")

cancer_ever_cox_ipcw <-
  coxph(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data_wide, weights = baseline_wt, cluster = id)

hr_1c <- cancer_ever_cox_ipcw %>% 
  tidy_hr() %>% 
  mutate(
    model = "IPTW + IPCW"
  )

rd_1c <- risks_km(cancer_ever_km_ipcw) %>% 
  mutate(
    model = "IPTW + IPCW"
  )

bind_rows(hr_1a, hr_1b, hr_1c)

bind_rows(rd_1a, rd_1b, rd_1c)

# 2. Time-varying cancer -----------------------------------------------------


# 2.1. IPTW for t-v cancer, with baseline covariates -------------------

smoke1 <- data_wide %>% select(id, smoke1, ht1, sbp1, bmi1, diabetes_prev)

data_long <- data_long %>% 
  left_join(smoke1)

cancer_den <-
  glm(
    cancer_v ~ bs(age_0, 3) + sex + education + apoe4 + as.factor(smoke1) + ht1 + 
      bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) + cohort,
    data = data_long,
    family = binomial
  )

summary(cancer_den)

cancer_num <- glm(cancer_v ~ 1, data = data_long)

summary(cancer_num)

data_long <- data_long %>%
  mutate(
    p_num = predict(cancer_num, type = "response"),
    p_denom = predict(cancer_den, type = "response"))
  
data_long <- data_long %>%
  group_by(id) %>% 
  mutate(w_cancer = ifelse(cancer_v == 1, p_num/p_denom, (1 - p_num)/(1- p_denom))) %>% 
  ungroup()

data_long %<>% 
  mutate(w_cancer_t = ifelse((w_cancer > quantile(w_cancer, 0.99)), 
                            quantile(w_cancer, 0.99), w_cancer))


# 2.2. IPCW. Weights on death - dementia ---------------------------------------------

death_den <- glm(
  competing_plr ~ cancer_v + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
    as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) +
    as.factor(smoke) + bs(sbp, 3) + bs(bmi, 3) + ht + ht_drug + hd_v + stroke_v + diab_v + cohort,
  data = data_long,
  family = quasibinomial
)

summary(death_den)

death_den %>% broom::tidy(exponentiate =TRUE) %>% View()

death_num <-
  glm(
    competing_plr ~ cancer_v + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
      as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) + cohort,
    data = data_long,
    family = binomial
  )

death_num %>% broom::tidy(exponentiate =TRUE) %>% View()

data_long$p_denom = predict(death_den, data_long, type = "response")

data_long$p_num = predict(death_num, data_long, type = "response")

data_long %<>%
  group_by(id) %>%
  mutate(w_death = cumprod(1 - p_num) / cumprod(1 - p_denom)) %>% 
  ungroup()

data_long %<>% 
  mutate(w_death_t = ifelse((w_death > quantile(w_death, 0.99)), 
                            quantile(w_death, 0.99), w_death))

data_long %>% ggplot(aes(x = as_factor(time), y = w_death_t)) +
  geom_boxplot()

data_long %<>%
  mutate(weights_both = w_cancer * w_death, 
         weights_both_t = w_cancer_t * w_death_t)
  
# 2.3. KM independent censoring, with confounding ------------------------------

km_unconditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = w_cancer_t)

plot_km(km_unconditional, "Net risk") + 
  labs(subtitle = "Unconditional independent censoring of death")


net_unconditional_cox <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = w_cancer_t)

tidy(net_unconditional_cox, exponentiate = TRUE)

risks_km(km_unconditional)


# 2.4. Independent censoring conditional on tv covariates ---------------

km_conditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = weights_both_t)

cox_conditional <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = weights_both_t)

tidy(cox_conditional, exponentiate = TRUE)

plot_km(km_conditional, "Direct effect") + 
  labs(subtitle = "With IPTW, conditional censoring on time-v")

risks_km(km_conditional)


# 3. Time to cancer -------------------------------------------------------


# 3.1. Weights for t2cancer -----------------------------------------------

# We want the probability of having a cancer dx at each 
# time point for each participant who has not yet being diagnosed

data_long %<>% 
  group_by(id) %>%
  mutate(pastcancer = ifelse(cancer_v == 1 & lag(cancer_v) == 0, 0, cancer_v),
         pastcancer = ifelse(is.na(pastcancer), 0, pastcancer)) %>% 
  ungroup() 


# denominator
mod <- glm(cancer_v ~ bs(age_0) + sex + education + apoe4 + as.factor(smoke1) + ht1 + 
             bs(sbp1, 3) + bs(bmi1, 3) + as.factor(diabetes_prev) + cohort,
           family = quasibinomial(), data = subset(data_long , 
                                                   pastcancer == 0))

tidy(mod, exponentiate = TRUE) %>% View()

data_long %<>%
  mutate(pred_a_den = ifelse(pastcancer == 1, 
                             1, # the pr of ever having cancer is 1 for those with cancer
                             predict(mod, type = 'response'))) 

data_long %>% ggplot(aes(pred_a_den)) + geom_histogram()

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
  mutate(sw_t2cancer = numcum/dencum,
         w_t2_cancer = 1/dencum)


data_long %<>% 
  mutate(sw_t2cancer_t = ifelse((sw_t2cancer > quantile(sw_t2cancer, 0.99)), 
                            quantile(sw_t2cancer, 0.99), sw_t2cancer))

# Weights
data_long %>% ggplot(aes(sw_t2cancer)) + geom_histogram()

# Multiply weights

data_long %<>%
  mutate(
    both_weights2 = w_death_t * sw_t2cancer_t)

data_long %>% ggplot(aes(x = as_factor(time), y = w_death_t)) +
  geom_boxplot()

# 3.2. KM independent censoring, with confounding ------------------------------

km_unconditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_t2cancer_t)

plot_km(km_unconditional, "Net risk") + 
  labs(subtitle = "Unconditional independent censoring of death")


net_unconditional_cox <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw_t2cancer_t)

tidy(net_unconditional_cox, exponentiate = TRUE)

risks_km(km_unconditional)


# 2.4. Independent censoring conditional on tv covariates ---------------

km_conditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = both_weights2)

cox_conditional <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = both_weights2)

tidy(cox_conditional, exponentiate = TRUE)

plot_km(km_conditional, "Direct effect") + 
  labs(subtitle = "With IPTW, conditional censoring on time-v")

risks_km(km_conditional)


# Bounds ------------------------------------------------------------------



