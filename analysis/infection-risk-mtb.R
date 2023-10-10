#### Libraries ####

library(tidyverse)
library(truncnorm)
library(LaplacesDemon)
source("utils/settings.R")
source("utils/plotting.R")

#### Data ###

# quanta
source("utils/quanta.R")

# CO2
co2 <- readRDS("data-clean/co2.rds")

# incidence
prev <- readRDS("data-clean/estimated-prevalence-mtb.rds")


#### Main analysis ####

main <- tibble(
  scenario = c("Low", "Medium", "High"),
  q = list(q_mtb_low, q_mtb_med, q_mtb_high)
) %>%
  mutate(scenario = factor(scenario, levels = c("Low", "Medium", "High"))) %>%
  unnest() %>%
  group_by(scenario) %>%
  mutate(draw = 1:n()) %>%
  left_join(prev %>%
              filter(group == "young population") %>%
              unnest() %>%
              group_by(country) %>%
              mutate(draw = 1:n()) %>%
              ungroup()) %>%
  left_join(data.frame(country = names(n.class), n = n.class)) %>%
  left_join(co2 %>% dplyr::select(country, f.mean)) %>%
  mutate(t = year) %>% 
  mutate(P = 1 - exp(-f.mean*q*I*t/n))

main_pl <- plot_risk(main) +
  ylab("Annual risk of infection (%, square-root scale)")

main_pl

ggsave("results/infection-risk-mtb-main.png", 
       width = 16, height = 10, units = "cm")
ggsave("results/infection-risk-mtb-main.eps", 
       width = 16, height = 10, units = "cm")


main %>%
  group_by(country, scenario) %>%
  summarize_at("P", list(q50 = median, q25 = q25, q75 = q75, q2.5 = q2.5, q97.5 = q97.5)) %>%
  mutate_if(is.numeric, function(x) round_k(x * 100, 1))

#### Outbreak ####

outbreak <- tibble(
  scenario = c("Low", "Medium", "High"),
  q = list(q_mtb_low, q_mtb_med, q_mtb_high)
) %>%
  mutate(scenario = factor(scenario, levels = c("Low", "Medium", "High"))) %>%
  unnest() %>%
  group_by(scenario) %>%
  mutate(draw = 1:n()) %>%
  ungroup() %>%
  left_join(prev %>%
              filter(group == "outbreak") %>%
              unnest() %>%
              group_by(country) %>%
              mutate(draw = 1:n()) %>%
              ungroup()) %>%
  left_join(data.frame(country = names(n.class), n = n.class)) %>%
  left_join(co2 %>% 
              dplyr::select(country, C.daily) %>% 
              unnest() %>% 
              dplyr::select(country, f.mean) %>%
              group_by(country) %>%
              sample_n(n.sample, T) %>%
              ungroup() %>%
              group_by(country) %>%
              mutate(draw = 1:n()) %>%
              ungroup()) %>%
  mutate(t = week) %>% 
  mutate(P = 1 - exp(-f.mean*q*I*t/n))

outbreak_pl <- plot_risk(outbreak) +
  ylab("Weekly risk of infection (%, square-root scale)")

outbreak_pl

ggsave("results/infection-risk-mtb-outbreak.png", 
       width = 16, height = 10, units = "cm") 
ggsave("results/infection-risk-mtb-outbreak.eps", 
      width = 16, height = 10, units = "cm") 

outbreak %>%
  group_by(country, scenario) %>%
  summarize_at("P", list(q50 = median, q25 = q25, q75 = q75, q2.5 = q2.5, q97.5 = q97.5)) %>%
  mutate_if(is.numeric, function(x) round_k(x * 100, 1))


#### Sens: Reported ####

sens.reported <- rbind(
  prev %>%
    filter(group == "young population") %>%
    unnest() %>%
    group_by(country) %>%
    mutate(draw = 1:n()) %>%
    ungroup(),
  prev %>%
    filter(group == "general population") %>%
    unnest() %>%
    group_by(country) %>%
    sample_n(n.sample, T) %>%
    mutate(draw = 1:n()) %>%
    ungroup()) %>%
  left_join(data.frame(q = q_mtb_med) %>%
              mutate(draw = 1:n())) %>%
  left_join(data.frame(country = names(n.class), n = n.class)) %>%
  left_join(co2 %>% dplyr::select(country, f.mean)) %>%
  mutate(t = year) %>% 
  mutate(P = 1 - exp(-f.mean*q*I*t/n)) %>%
  rename(scenario = group) %>%
  mutate(scenario = factor(scenario, levels = c("general population", "young population")))

sens.reported_pl <- plot_risk(sens.reported) +
  scale_x_discrete(labels = c("General population", "Young population")) +
  xlab("Prevalence") +
  ylab("Annual risk of infection (%, square-root scale)")

sens.reported_pl

ggsave("results/infection-risk-mtb-sensitivity_population.png", 
       width = 12, height = 10, units = "cm")
ggsave("results/infection-risk-mtb-sensitivity_population.eps", 
       width = 12, height = 10, units = "cm")


#### Sens: Outdoor CO2 ####

sens.Co <- tibble(
  scenario = c("Low", "Medium", "High"),
  q = list(q_mtb_low, q_mtb_med, q_mtb_high)
) %>%
  mutate(scenario = factor(scenario, levels = c("Low", "Medium", "High"))) %>%
  unnest() %>%
  group_by(scenario) %>%
  mutate(draw = 1:n()) %>%
  left_join(prev %>%
              filter(group == "young population") %>%
              unnest() %>%
              group_by(country) %>%
              mutate(draw = 1:n()) %>%
              ungroup()) %>%
  left_join(data.frame(country = names(n.class), n = n.class)) %>%
  left_join(co2 %>% dplyr::select(country, f.mean.sens)) %>%
  mutate(t = year) %>% 
  mutate(P = 1 - exp(-f.mean.sens*q*I*t/n))

sens.Co <- sens.Co %>%
  filter(scenario == "Medium") %>%
  mutate(scenario = country,
         type = "600ppm") %>% 
  rbind(rbind(
    main %>%
      filter(scenario == "Medium") %>%
      filter(country != "South Africa") %>%
      mutate(scenario = "South Africa"),
    main %>% 
      filter(scenario == "Medium") %>%
      filter(country != "Switzerland") %>%
      mutate(scenario = "Switzerland"),
    main %>% 
      filter(scenario == "Medium") %>%
      filter(country != "Tanzania") %>%
      mutate(scenario = "Tanzania")) %>% 
      mutate(type = "400ppm"))


sens.Co_pl <- plot_risk.sens(sens.Co) +
  ylab("Annual risk of infection (%, square-root scale)")

sens.Co_pl

ggsave("results/infection-risk-mtb-sensitivity_outdoorCO2.png", 
       width = 16, height = 11, units = "cm")  
ggsave("results/infection-risk-mtb-sensitivity_outdoorCO2.eps", 
       width = 16, height = 11, units = "cm")  
