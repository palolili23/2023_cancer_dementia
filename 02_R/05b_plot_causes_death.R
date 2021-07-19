library(rio)
library(tidyverse)

ggthemes::scale_fill_colorblind()

wide_data <- import(here::here("01_data", "clean_data", "wide_after_truncation.RData"))

cases_death <- import(here::here("01_data", "raw_data", "causes_death.RData"))

data <- wide_data %>% 
  left_join(cases_death, by = c("id" = "ergoid"))

# data %>%
#   group_by(cancer) %>%
#   count(dementia) %>% 
# mutate(prop = paste0(round(100*(n/sum(n)),0), "%")) %>% 
#   gt::gt()

count_cancer_death <- data %>% 
  filter(competing_plr == 1) %>% 
  group_by(cancer_v) %>% 
  count(category) %>% 
  mutate(prop = 100*(n/sum(n))) 

causes_death <- count_cancer_death %>% 
  filter(category != "Alive") %>% 
  mutate(category = str_to_title(category)) %>% 
  mutate(cancer = ifelse(cancer_v == 0, "No cancer", "Incident cancer"),
    category = as_factor(category),
         category = fct_reorder(category, prop),
         cancer = as_factor(cancer)) %>% 
  ggplot(aes(cancer, prop, fill = category)) +
  geom_col() +
  ggthemes::scale_fill_tableau() +
  coord_flip() + 
  labs(y = "% of participants that died prior to a dementia diagnosis",
       x = NULL,
       fill = NULL) +
  theme_minimal(base_family = "serif") +
  theme(legend.position = "bottom",
        legend.text = element_text(size=12)) +
  theme(strip.text.x = element_text(size = 11),
        strip.background = element_rect(fill=NA),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggsave(filename = "causes_death.tiff",
       plot = causes_death,
       path = here::here("03_figs"),
       device = "tiff",
       width = 8,
       height = 4.1,
       dpi = "retina")


count_dementia_death <- data %>% 
  filter(competing_plr == 0) %>% 
  group_by(outcome_plr) %>% 
  count(category) %>% 
  mutate(prop = 100*(n/sum(n))) 

causes_death_dementia <- count_dementia_death %>% 
  # filter(category != "Alive") %>% 
  mutate(category = str_to_title(category)) %>% 
  mutate(outcome_plr = ifelse(outcome_plr == 0, "Had no dementia diagnosis", "Had dementia diagnosis"),
         category = as_factor(category),
         # category = fct_reorder(category, prop),
         outcome_plr = as_factor(outcome_plr)) %>% 
  ggplot(aes(outcome_plr, prop, fill = category)) +
  geom_col() +
  ggthemes::scale_fill_tableau() +
  coord_flip() + 
  theme_minimal() +
  labs(y = "% of participants that died",
       x = NULL,
       fill = NULL) +
  theme(legend.position = "bottom",
        legend.text = element_text(size=10)) +
  theme(
    strip.text.x = element_text(size = 10),
    # strip.background = element_blank(),
    strip.background = element_rect(fill=NA),
    axis.text=element_text(size=10),
    axis.title=element_text(size=10))
