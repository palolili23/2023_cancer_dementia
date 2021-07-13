library(tidyverse)
ever <- rio::import(here::here("02_R", "rd_results_ever.csv"))
tv <- rio::import(here::here("02_R", "rd_results_tv.csv"))

hr <-rio::import(here::here("02_R", "hr_results.csv")) %>%
  mutate(Proxy = ifelse(Proxy == "Time to cancer", "Time-varying", Proxy)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(Proxy = str_replace(Proxy, "Ever vs. Never", "Ever vs. never")) %>% 
  mutate(`Hazard Ratio` = paste0(point, " (", conf.low, ",", conf.high, ")")) %>% 
  select(Proxy, model, `Hazard Ratio`)

rr <- ever %>%
  bind_rows(tv) %>%
  select(Proxy, model, rr, conf) %>%
  pivot_wider(names_from = conf,
              values_from = rr) %>%
  mutate(contrast = "Risk Ratio") %>%
  rename(conf.high = conf.higher) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(`Risk Ratio` = paste0(point, " (", conf.low, ",", conf.high, ")")) %>% 
  select(Proxy, model, `Risk Ratio`)


rd <- ever %>%
  bind_rows(tv) %>%
  select(Proxy, model, rd, conf) %>%
  pivot_wider(names_from = conf,
              values_from = rd) %>%
  mutate(contrast = "Risk Difference") %>%
  rename(conf.high = conf.higher) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  mutate(`Risk Difference`= paste0(point, " (", conf.low, ",", conf.high, ")")) %>% 
  select(Proxy, model, `Risk Difference`)

complete_table <- rd %>% left_join(rr) %>% left_join(hr)

complete_table <- complete_table %>% 
  mutate(Proxy = ifelse(Proxy == "Ever vs. never", "A", "B")) %>%
  rename(Scenario = Proxy) %>% 
  mutate(model = 
           case_when(model == "IPTW" ~ "With weights for confounding",
                     model == "IPTW + IPCW" ~ "With weights for confounding and for censoring for death",
                     TRUE ~ model),
         model = as_factor(model),
         model = fct_relevel(model, c("Unadjusted", "With weights for confounding", "With weights for confounding and for censoring for death"))) %>% 
  group_by(Scenario) %>% 
  arrange(model)

bounds <- rio::import(here::here("02_R", "table_results_bounds.csv"))

bounds_tv_ipw <- bounds %>% slice(4,6) %>% select(model, rr) %>% pull(rr)

bounds_tv_ipw <- paste0("RR: ", bounds_tv_ipw[1], ", ", bounds_tv_ipw[2])


