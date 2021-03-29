library(rio)
library(tidyverse)
library(wesanderson)
ggthemes::scale_fill_colorblind()

wide_data <- import(here::here("01_data", "clean_data", "wide_noltfu.RData"))

cases_death <- import(here::here("01_data", "raw_data", "causes_death.RData"))

data <- wide_data %>% 
  left_join(cases_death, by = c("id" = "ergoid"))

data %>%
  group_by(cancer) %>%
  count(dementia) %>% 
mutate(prop = paste0(round(100*(n/sum(n)),0), "%")) %>% 
  gt::gt()

count_cancer_death <- data %>% 
  filter(cohort == 1, dementia == 2) %>%
  group_by(cancer) %>% 
  count(category) %>% 
  mutate(prop = 100*(n/sum(n))) 

count_cancer_death %>% 
  filter(category != "Alive") %>% 
  mutate(cancer = ifelse(cancer == 0, "No cancer", "Incident cancer"),
    category = as_factor(category),
         category = fct_reorder(category, prop),
         cancer = as_factor(cancer)) %>% 
  ggplot(aes(cancer, prop, fill = category)) +
  geom_col() +
  ggthemes::scale_fill_tableau() +
  coord_flip() + 
  theme_minimal() +
  labs(y = "% of participants that died",
       x = NULL,
       fill = NULL)
    
