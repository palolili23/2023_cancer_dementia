library(ggplot2)
library(tidyverse)
library(rio)
library(wesanderson)
library(survival)
library(survminer)
library(broom)
library(splines)
library(magrittr)

# Import data -------------------------------------------------------------
source(here::here("02_R", "04b_auxiliary_fx.R"))

data_long <-
  import(here::here("01_data", "clean_data", "dementia_long.RData"))

data_wide <-
  data_long %>%
  group_by(id) %>% 
  slice(n()) %>% 
  ungroup()

## Use these numbers for chart

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

# Weights for death -------------------------------------------------------

death_den <- glm(competing_plr ~ 
                   cancer_v + bs(time,3) + bs(age_0,3) + sex + education + apoe4 + 
                   bs(sbp,3) + bs(bmi,3 ) + hd_v + stroke_v + diab_v, 
                 data = data_long, family = quasibinomial)

summary(death_den)

death_den %>% broom::tidy() %>% View()

death_num <- glm(competing_plr ~ 1, data = data_long, family = binomial)

data_long$p_denom = predict(death_den, data_long, type = "response")

data_long$p_num = predict(death_num, data_long, type = "response")

data_long %<>%
  group_by(id) %>%
  mutate(sw = cumprod(1 - p_num) / cumprod(1 - p_denom),
  ) %>% 
  ungroup()

data_long %<>% 
  mutate(sw = ifelse((
    sw > quantile(sw, 0.99)
  ), quantile(sw, 0.99), sw))

data_long %>% ggplot(aes(x = as_factor(time), y = sw)) +
  geom_boxplot()


# 1. EVER VS. NEVER -------------------------------------------------------

# 1.1. Crude risks --------------------------------------------------------

cancer_ever <-
  survfit(Surv(t2dem, as.factor(outcome_plr)) ~ cancer_v, data_wide)

survminer::ggcompetingrisks(cancer_ever) +
  ggthemes::scale_fill_tableau() +
  theme_minimal() +
  theme(
    legend.position = "bottom")

plot_cif(cancer_ever, "Crude risk")

## subdistribution cox

data_wide %<>% 
  mutate(t2dem_no_death = 
           ifelse(dementia_20 == 2, 240, t2dem_20))

cancer_ever_subcox <-
  coxph(Surv(t2dem_no_death, dementia_20 == 1) ~ cancer_20, data_wide)

## risks and hazards
tidy(cancer_ever_subcox, exponentiate = TRUE)

risks_cif(cancer_ever)

# 1.2. Net risks ----------------------------------------------------------

# 1.2.a. Independent censoring --------------------------------------------

cancer_ever_km <-
  survfit(Surv(t2dem_20, dementia_20 == 1) ~ cancer_20, data_wide)

plot_km(cancer_ever_km, "Net risks") +
  labs(subtitle = "Under unconditional independent censoring of death")

cancer_ever_cox <-
  coxph(Surv(t2dem_20, dementia_20 == 1) ~ cancer_20, data_wide)

tidy(cancer_ever_cox, exponentiate = TRUE)

risks_km(cancer_ever_km)


# 1.2.b. Independent censoring conditional on tv covariates ---------------

cancer_ever_km_tv <-
  survfit(Surv(
    tstart, time2 = fuptime, event = outcome_plr) ~ cancer_20,
    data = data_long, cluster = id, weights = sw)

plot_km(cancer_ever_km_tv, "Net risks") +
  labs(subtitle = "Independent censoring of death conditional on tv covs")

cancer_ever_cox_tv <-
  coxph(Surv(
    tstart, time2 = fuptime, event = outcome_plr) ~ cancer_20,
    data = data_long, cluster = id, weights = sw)

tidy(cancer_ever_cox_tv, exponentiate = TRUE)

risks_km(cancer_ever_km_tv)


# 2. PREVALENT AND INCIDENT CANCER -------------------------------------------

# 2.1. Crude risk under time-varying cancer diagnosis ---------------------

cancer_incident_cif <- survfit(Surv(tstart, time2 = fuptime, 
                    event = as_factor(both_outcomes)) ~ cancer_v,
               data = data_long, id = id)

plot_cif(cancer_incident_cif, "Crude risk, time-varying cancer")

cancer_incident_cox <-
  coxph(
    Surv(tstart, time2 = fuptime, event = as_factor(both_outcomes)) ~ cancer_v,
    data = data_long,
    id = id
  )

tidy(cancer_incident_cox, exponentiate = TRUE)


# 2.2. Net risks with incident cancer -------------------------------------


# 2.2.a. Independent censoring --------------------------------------------

km_unconditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
    data = data_long, cluster = id)

plot_km(km_unconditional, "Net risk") + 
  labs(subtitle = "Unconditional independent censoring of death")


net_unconditional_cox <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id)

tidy(net_unconditional_cox, exponentiate = TRUE)

risks_km(km_unconditional)


# 1.2.b. Independent censoring conditional on tv covariates ---------------

km_conditional <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw)

cox_conditional <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = sw)

tidy(cox_conditional, exponentiate = TRUE)

plot_km(km_conditional, "Direct effect") + 
  labs(subtitle = "Without IPTW, conditional censoring on time-v")

risks_km(km_conditional)


#3. INCIDENT CANCER WITH IPW ------------------------------------------------


# 3.1. Weights calculation ------------------------------------------------

data_long %<>% 
  group_by(id) %>%
  mutate(pastcancer = ifelse(cancer_v == 1 & lag(cancer_v) == 0, 0, cancer_v),
         pastcancer = ifelse(is.na(pastcancer), 0, pastcancer)) %>% 
  ungroup() 

## Fitting weights for cancer
mod <- glm(cancer_v ~ bs(time, 3) + sex + age_0 +
             cohort + education + smoke + bs(sbp, 3) +
            bmi + I(bmi^2) + stroke_v + diab_v + hd_v,
           family = quasibinomial(), data = subset(data_long , 
                                              pastcancer == 0))

tidy(mod, exponentiate = TRUE) %>% gt::gt()

data_long %<>%
  mutate(pred_a_den = 
           predict(mod, new = data_long, type = "response"))

mod_num <- glm(cancer_v ~ bs(time, 3) + sex + age_0 +
                 cohort + education,
               family = quasibinomial(),
               data = subset(data_long, pastcancer == 0))

data_long %<>% 
  mutate(pred_a_num = 
           predict(mod_num, new = data_long, type = "response"))

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
    both_weights = sw_cancer*sw)


# 3.2. Crude risks with weights for cancer -------------------------------------
cancer_incident_cif_ipw <- survfit(
  Surv(tstart, time2 = fuptime,
       event = as_factor(both_outcomes)) ~ cancer_v,
  data = data_long,
  id = id,
  weights = sw_cancer
)

plot_cif(cancer_incident_cif_ipw,
         "Crude risk, IPW time-varying cancer")

cancer_incident_cox_ipw <-
  coxph(
    Surv(tstart, time2 = fuptime, event = as_factor(both_outcomes)) ~ cancer_v,
    data = data_long,
    id = id,
    weights = sw_cancer
  )

tidy(cancer_incident_cox_ipw, exponentiate = TRUE)

# 3.3. Net risk with weights for cancer and death ------------------------------------

km_ipweighted <- survfit(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = both_weights)

cox_ipweighted <- coxph(Surv(
  tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
  data = data_long, cluster = id, weights = both_weights)

tidy(cox_ipweighted, exponentiate = TRUE)

plot_km(km_ipweighted, "Direct effect") + 
  labs(subtitle = "With IPTW, conditional censoring on time-v")

risks_km(km_ipweighted)

