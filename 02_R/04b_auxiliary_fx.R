  mytheme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    plot.caption = element_text(hjust = 0))



risks_cif <- function(model){
  model %>% 
    broom::tidy() %>%
    filter(state == "1") %>%
    group_by(strata) %>% 
    slice(n()) %>% 
    select(estimate, strata) %>% 
    pivot_wider(names_from = strata, values_from = estimate) %>% 
    mutate(rd = .[[2]] - .[[1]],
           rr = .[[2]] / .[[1]]) %>% 
    mutate_at(c(1:3), ~.*100)
}



plot_cif <- function(model, title, ...){
  
  label_exp <- model %>% 
    tidy() %>% slice(1) %>% pull(strata)
  
  tidy <- model %>% 
    broom::tidy() %>% 
    filter(state == "1") %>%
    select(time, strata, estimate, conf.high, conf.low) %>%
    rename(CIF = estimate) %>% 
    mutate(strata = ifelse(strata == label_exp, "Free of cancer", "Incident cancer"))
  
  plot <- tidy %>% 
    ggplot(aes(time, CIF, group = strata)) +
    geom_line(aes(color = strata), size = 0.7) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high) ,alpha = 0.2) +
    scale_color_manual(values = c("#011A5E", "#e4a803")) +
    scale_y_continuous(limits = c(0, 0.70)) +
    labs(
      title = paste0(title),
      color = NULL,
      y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(plot)
  
}


# Km ----------------------------------------------------------------------

plot_km <- function(model, title, ...){
  
  label_exp <- model %>% 
    tidy() %>% slice(1) %>% pull(strata)
  
  tidy <- model %>% 
    broom::tidy() %>% 
    transmute(
      time = time, 
      strata = ifelse(strata == label_exp, 
                               "Free of cancer", "Incident cancer"),
      cif = 1 - estimate,
      conf.low2 = 1- conf.low,
      conf.high2 = 1 - conf.high) %>% 
    rename(conf.high = conf.low2,
           conf.low = conf.high2)
  
  plot <- tidy %>% 
    ggplot(aes(time, cif, group = strata)) +
    geom_line(aes(color = strata), size = 0.7) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high) ,alpha = 0.2) +
    scale_color_manual(values = c("#011A5E", "#e4a803")) +
    scale_y_continuous(limits = c(0, 0.70)) +
    labs(
      title = paste0(title),
      color = NULL,
      y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(plot)
  
}

risks_km <- function(model){
  model %>% 
    broom::tidy() %>%
    group_by(strata) %>% 
    slice(n()) %>% 
    mutate(estimate = 1-estimate) %>% 
    select(estimate, strata) %>% 
    pivot_wider(names_from = strata, values_from = estimate) %>% 
    mutate(rd = .[[2]] - .[[1]],
           rr = .[[2]] / .[[1]]) %>% 
    mutate_at(c(1:3), ~.*100)} 

tidy_hr <- function(model) {
 model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, contains("conf"))
}
