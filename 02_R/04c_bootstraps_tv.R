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

## If competing_plr = 1, then outcome_plr = NA

data_long %<>% 
  mutate(outcome_plr = ifelse(competing_plr == 1, NA, outcome_plr))

## Merge with baseline
smoke1 <- data_wide %>% select(id, smoke1, ht1, sbp1, bmi1)

data_long <- data_long %>% 
  left_join(smoke1)

# Weighted KM function for bootstrap  -------------------------------------------------------


km_tv <- function(data_long, crude = TRUE, ipcw = FALSE) {
  if (crude != FALSE) {
    
    model <- survfit(Surv(
        tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
        data = data_long)
    
    
    output <- risks_km(model) %>% 
      mutate(
        model = "Unadjusted")
  }
  
  else {  
    
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
    
    
    
    if (ipcw != FALSE) {
      death_den <- glm(
        competing_plr ~ cancer_lag + bs(time, 3) + bs(age_0, 3) + sex + education + apoe4 +
          as.factor(smoke_lag) + bs(sbp_lag, 3) + bs(bmi_lag, 3) + 
          ht_lag + ht_drug_lag + hd_lag + stroke_lag + diab_lag + cohort,
        data = data_long,
        family = quasibinomial
      )
      
      summary(death_den)
      
      death_num <-
        glm(
          competing_plr ~ cancer_lag + bs(time, 3) + cohort,
          data = data_long,
          family = binomial
        )
      
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
      
      data_long %<>%
        mutate(weights_both = w_cancer * w_death, 
               weights_both_t = w_cancer_t * w_death_t,
               sw_weights_both = sw_cancer * sw_death,
               sw_weights_both_t = sw_cancer_t * sw_death_t)
      
      model <- survfit(Surv(
        tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
        data = data_long, cluster = id, weights = sw_weights_both_t)
      
      output <- risks_km(model) %>%
        mutate(model = "IPTW + IPCW")
    }
    else{

      model <- survfit(Surv(
        tstart, time2 = fuptime, event = outcome_plr) ~ cancer_v,
        data = data_long, cluster = id, weights = sw_cancer_t)
      
      output <- risks_km(model) %>% 
        mutate(model = "IPTW")
    }
    
  }
  
  return(output)  
}


# Bootstrap function ------------------------------------------------------

risks_boot_long <- function(data_long, n, seed, crude = TRUE, ipcw = FALSE){
  # Set seed
  set.seed(seed)
  id_vector <-data_long %>% distinct(id) %>% pull(id)
  n_size <- length(id_vector)
  
  # Creates bootsamples and runs the model to each sample
  bootsamps <- replicate(n = n, expr = {
    d <- sample(id_vector,size = n_size, replace = T)
    ds_b <- data_long %>% filter(id %in% d)
    output <- km_tv(ds_b, crude, ipcw) 
    return(output)
  }, simplify = F)
  
  totalboot <- bind_rows(bootsamps)
  
  totalboot <- totalboot %>% 
    summarise(across(.cols = c(1:4), 
                     ~quantile(.x, probs = c(0.025, 0.975))))
  
  return(totalboot)
}

# Results -----------------------------------------------------------------

point_km_crude <- km_tv(data_long, crude = TRUE)
point_km_crude$conf <- "point"

point_km_iptw <- km_tv(data_long, crude = FALSE, ipcw = FALSE)
point_km_iptw$conf <- "point"

point_km_ipcw <- km_tv(data_long, crude = FALSE, ipcw = TRUE)
point_km_ipcw$conf <- "point"

boots_km_crude <- risks_boot_long(data_long, seed = 123, n = 100, crude = TRUE, ipcw = FALSE)

boots_km_crude$conf <- c("conf.low", "conf.higher")
boots_km_crude$model <- "Unadjusted"

boots_km_iptw <- risks_boot_long(data_long, seed = 123, n = 100, crude = FALSE, ipcw = FALSE)

boots_km_iptw$conf <- c("conf.low", "conf.higher")
boots_km_iptw$model <- "IPTW"

boots_km_ipcw <- risks_boot_long(data_long, seed = 123, n = 100, crude = FALSE, ipcw = TRUE)
boots_km_ipcw$conf <- c("conf.low", "conf.higher")
boots_km_ipcw$model <- "IPTW + IPCW"


rd_results_tv <- bind_rows(point_km_crude, point_km_iptw, point_km_ipcw) %>%
  bind_rows(boots_km_crude, boots_km_iptw, boots_km_ipcw) %>%
  mutate(Proxy = "Time-varying")

export(rd_results_tv, here::here("02_R", "rd_results_tv.csv"))
