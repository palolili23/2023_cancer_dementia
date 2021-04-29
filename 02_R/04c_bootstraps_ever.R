source(here::here("02_R", "04b_auxiliary_fx.R"))

data_wide <-
  import(here::here("01_data", "clean_data", "wide_after_truncation.RData"))

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

# Weighted KM function for bootstrap ------------------------------------------------------

km_ever <- function(data, crude = TRUE, ipcw = FALSE) {
  if (crude != FALSE) {
    
    model <-
      survfit(Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v, data)
    
    output <- risks_km(model) %>% 
      mutate(
        model = "Unadjusted"
      )
  }
  
  else {  
    cancer_den <-
      glm(
        cancer_v ~ bs(age_0) + sex + education + apoe4 + as.factor(smoke1),
        data = data,
        family = binomial
      )
    
    cancer_num <- glm(cancer_v ~ 1, data = data)
    
    data <- data %>%
      mutate(
        p_num = predict(cancer_num, type = "response"),
        p_denom = predict(cancer_den, type = "response"),
        w_cancer = ifelse(cancer_v == 1, p_num/p_denom, (1 - p_num)/(1- p_denom)))
    
    data %<>%
      mutate(w_cancer_t = ifelse(
        w_cancer > quantile(w_cancer, 0.99),
        quantile(w_cancer, 0.99),
        w_cancer
      ))
    
    if (ipcw != FALSE) {
      death_den <-
        glm(
          competing_plr ~ bs(age_0, 3) + sex + education + apoe4 + 
            as.factor(smoke1) + ht1 + bs(sbp1, 3) + bs(bmi1,3) +
            as.factor(diabetes_prev) +  cancer_v + cohort,
          data = data,
          family = binomial
        )
      
      death_num <- glm(competing_plr ~ 1, data = data)
      
      data <- data %>%
        mutate(
          dp_num = predict(death_num, type = "response"),
          dp_denom = predict(death_den, type = "response"),
          w_death = if_else(competing_plr == 0, (1-dp_num)/(1-dp_denom), 1))
      
      data %<>%
        mutate(w_death_t = ifelse(
          w_death > quantile(w_death, 0.99),
          quantile(w_death, 0.99),
          w_death
        ))
      
      data %<>%
        mutate(
          baseline_w = w_cancer * w_death,
          baseline_wt = w_cancer_t * w_death_t)
      
      
      model <-
        survfit(
          Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
          data,
          weights = baseline_wt,
          cluster = id
        )
      
      output <- risks_km(model) %>%
        mutate(model = "IPTW + IPCW")
    }
    else{
      model <-
        survfit(
          Surv(t2dem_efu, dementia_efu == 1) ~ cancer_v,
          data,
          weights = w_cancer_t,
          cluster = id
        )
      
      output <- risks_km(model) %>% 
        mutate(model = "IPTW")
    }
    
  }
  
  return(output)  
}


# Bootstrap function ------------------------------------------------------

risks_boots <- function(data, n,
                        seed = 123,
                        crude = TRUE, ipcw = FALSE){
  # Set seed
  set.seed(seed)
  
  # Creates bootsamples and runs the model to each sample
  bootsamps <- replicate(n = n, expr = {
    d <- sample(1:nrow(data),size = nrow(data), replace = T)
    ds_b <- data[d,] %>%
      mutate(id = row_number()) 
    output <- km_ever(ds_b, crude, ipcw)
    return(output)
  }, simplify = F)
  
  totalboot <- bind_rows(bootsamps)
  
  totalboot <- totalboot %>% 
    summarise(across(.cols = c(1:4), 
                     ~quantile(.x, probs = c(0.025, 0.975))))
  
  return(totalboot)
}

# Bootstrap ---------------------------------------------------------------

point_estimate_crude <- km_ever(data_wide, crude = TRUE, ipcw = TRUE)
point_estimate_crude$conf <- "point"

point_estimate_iptw <- km_ever(data_wide, crude = FALSE, ipcw = FALSE)
point_estimate_iptw$conf <- "point"

point_estimate_ipcw <- km_ever(data_wide, crude = FALSE, ipcw = TRUE)
point_estimate_ipcw$conf <- "point"

km_crude_boots <- risks_boots(data_wide, n = 500, seed = 123, crude = TRUE, ipcw = FALSE)
km_crude_boots$conf <- c("conf.low", "conf.higher")
km_crude_boots$model <- "Unadjusted"

km_ipw_boots <- risks_boots(data_wide, n = 500, seed = 123, crude = FALSE, ipcw = FALSE)
km_ipw_boots$conf <- c("conf.low", "conf.higher")
km_ipw_boots$model <- "IPTW"

km_ipcw_boots <- risks_boots(data_wide, n = 500, seed = 123, crude = FALSE, ipcw = TRUE)
km_ipcw_boots$conf <- c("conf.low", "conf.higher")
km_ipcw_boots$model <- "IPTW + IPCW"

rd_results_ever <- bind_rows(point_estimate_crude, point_estimate_ipcw, point_estimate_iptw) %>%
  bind_rows(km_crude_boots, km_ipw_boots, km_ipcw_boots) %>%
  mutate(Proxy = "Ever vs. never")

export(rd_results_ever, here::here("02_R", "rd_results_ever.csv"))

