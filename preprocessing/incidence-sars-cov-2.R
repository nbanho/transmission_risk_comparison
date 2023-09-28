#### Libraries ####

library(tidyverse)
library(lubridate)
library(LaplacesDemon)
source("utils/settings.R")

#### Excess deaths ####

excess <- read.csv("data-raw/export_country_per_100k.csv")  %>% 
  filter(iso3c %in% c("ZAF", "CHE", "TZA")) %>% 
  mutate(date = ymd(date)) %>%        # Convert date column to date format
  filter(date <= ymd("2021-01-31") & date >= ymd("2020-03-01")) %>% 
  mutate(period = case_when(
    date >= ymd("2020-03-01") & date <= ymd("2020-05-31") ~ "April",
    date >= ymd("2020-06-01") & date <= ymd("2020-08-31") ~ "July",
    date >= ymd("2020-09-01") & date <= ymd("2020-11-30") ~ "October",
    date >= ymd("2020-12-01") & date <= ymd("2021-01-31") ~ "January",
    TRUE ~ NA_character_),
    country = case_when(iso3c == "CHE" ~ "Switzerland",
                        iso3c == "ZAF" ~ "South Africa",
                        iso3c == "TZA" ~ "Tanzania")) %>% 
  mutate(excess_deaths = estimated_daily_excess_deaths_per_100k * 7) %>% 
  dplyr::select(country, period, excess_deaths) %>% 
  nest(data = excess_deaths) %>%
  rename(excess_deaths = data) %>%
  mutate(excess_deaths = map(excess_deaths, unlist))

#### IFR ####

ifr <- tibble(country = c(rep("South Africa", 4), rep("Switzerland", 4), rep("Tanzania", 4)),
              period = rep(c("April", "July", "October", "January"), 3),
              IFR_median = c(0.511, 0.400, 0.348, 0.331, 1.687, 1.322, 1.164, 1.113, 0.169, 0.133, 0.116, 0.111),
              IFR_q2.5 = c(0.199, 0.173, 0.145, 0.137, 1.202, 0.943, 0.835, 0.766, 0.071, 0.060, 0.057, 0.052),
              IFR_q97.5 = c(1.275, 0.987, 0.774, 0.724, 2.782, 2.148, 1.743, 1.621, 0.339, 0.256, 0.218, 0.193)) %>% 
  mutate_at(vars(starts_with("IFR")), ~ . / 100) 


#### IFR-based simulation approach ####

# class size
nclass <- data.frame(country = ctry_ord, n = n.class)

# combine
incidence <- left_join(excess, ifr) %>% left_join(nclass)  %>%
  mutate(country = factor(country, levels = ctry_ord),
         period = factor(period, levels = c("April", "July", "October", "January")))

# function to simulate
sim_inci_ifr <- function(n_samples, exc_deaths, med_ifr, low_ifr, high_ifr) {
  
  # simulate excess deaths
  r_exc_deaths <- sample(exc_deaths, n_samples, T)
  
  # simulate ifr
  r_ifr <- rtrunc(n_samples, "norm", a = 0,
                  mean = med_ifr,
                  sd = (high_ifr - low_ifr) / (2 * qnorm(0.975)) )
  
  # compute incidence
  r_inci <- r_exc_deaths / r_ifr
  
  # set incidence to 0 if negative
  r_inci <- ifelse(r_inci < 0, 0, r_inci)
 
  return(r_inci)
}

incidence <- incidence %>%
  mutate(inci = pmap(list(n.sample, excess_deaths, IFR_median, IFR_q2.5, IFR_q97.5), sim_inci_ifr))


#### Summaries IFR incidence ####

# excess deaths
incidence %>%
  dplyr::select(country, period, matches("excess")) %>%
  unnest(excess_deaths) %>%
  group_by(country, period) %>%
  summarize(across(excess_deaths, list(median, q2.5, q97.5))) %>%
  ungroup() %>%
  arrange(country, period) %>%
  mutate_if(is.numeric, round_k, 2)
  
# IFR
incidence %>%
  dplyr::select(country, period, matches("IFR")) %>%
  arrange(country, period) %>%
  mutate_if(is.numeric, round_k, 3)

# incidence by period
incidence %>%
  dplyr::select(country, period, inci, n) %>%
  unnest(c(inci)) %>%
  unnest(c(inci)) %>%
  group_by(country, period) %>%
  summarize(across(inci, list(median, q2.5, q97.5))) %>%
  ungroup() %>%
  arrange(country, period) %>%
  mutate_if(is.numeric, round_k, 0)

# incidence overall
incidence %>%
  dplyr::select(country, period, inci, n) %>%
  unnest(c(inci)) %>%
  unnest(c(inci)) %>%
  group_by(country, period) %>%
  mutate(draw = 1:n()) %>%
  ungroup() %>%
  group_by(country, draw) %>%
  summarize(inci = mean(inci)) %>%
  ungroup() %>%
  group_by(country) %>%
  summarize_at("inci", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 0)
  
# I overall
incidence %>%
  dplyr::select(country, period, inci, n) %>%
  unnest(c(inci)) %>%
  unnest(c(inci)) %>%
  mutate(I = inci / 1e5 * n) %>%
  group_by(country, period) %>%
  mutate(draw = 1:n()) %>%
  ungroup() %>%
  group_by(country, draw) %>%
  summarize(I = mean(I)) %>%
  ungroup() %>%
  group_by(country) %>%
  summarize_at("I", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 3)


#### Reported incidence ####

# reported cases in South Africa in the general population
reported_sa <- read.csv("data-raw/WHO-COVID-19-global-data.csv") %>% 
  rename(date = Date_reported,
         country = Country) %>%
  mutate(date = ymd(date), # Convert Date_reported to a date object
         week = floor_date(date, unit = "week")) %>%  
  filter(country == "South Africa",
         date >= ymd("2020-03-01") & date <= ymd("2021-01-31")) %>%
  group_by(country, week) %>%
  summarize(inci = sum(New_cases) / 59390000 * 1e5) %>%
  ungroup() 
  

# reported cases in Switzerland in the general population
reported_ch <- read.csv("data-raw/COVID19CasesRawData_AKL10_d.csv") %>% 
  filter(date >= ymd("2020-03-01") & date <= ymd("2021-01-31")) %>% 
  group_by(date) %>% 
  summarise(daily_cases = sum(entries)) %>% 
  mutate(date = ymd(date),
         week = floor_date(date, unit = "week")) %>% 
  group_by(week) %>% 
  summarise(inci = sum(daily_cases) / 8700000 * 1e5) %>%
  ungroup() %>%
  mutate(country = "Switzerland") %>%
  dplyr::select(country, inci, week)

# reported cases in Switzerland in the young population
reported_ch_young <- read.csv("data-raw/COVID19CasesRawData_AKL10_d.csv") %>% 
  filter(ageRange == "10 - 19",
         date >= ymd("2020-03-01") & date <= ymd("2021-01-31")) %>% 
  group_by(date) %>% 
  summarise(daily_cases = sum(entries)) %>% 
  mutate(date = ymd(date),
         week = floor_date(date, unit = "week")) %>% 
  group_by(week) %>% 
  summarise(inci = sum(daily_cases) / 1755000 * 1e5) %>%
  ungroup() %>%
  mutate(country = "Switzerland") %>%
  dplyr::select(country, inci, week)


#### Summaries reported incidence ####

# reported incidence in the general population of South Africa
reported_sa %>%
  group_by(country) %>%
  summarize_at("inci", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 0)

# reported incidence in the general population of Switzerland
reported_ch %>%
  group_by(country) %>%
  summarize_at("inci", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 0)

# reported incidence in the young population of Switzerland
reported_ch_young %>%
  group_by(country) %>%
  summarize_at("inci", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 0)


#### Assumptions: I ####

# IFR-based approach
I_IFR <- incidence %>%
  dplyr::select(country, period, inci, n) %>%
  unnest(c(inci)) %>%
  unnest(c(inci)) %>%
  mutate(I = inci / 1e5 * n) %>%
  group_by(country, period) %>%
  mutate(draw = 1:n()) %>%
  ungroup() %>%
  group_by(country, draw) %>%
  summarize(I = mean(I)) %>%
  ungroup() %>%
  dplyr::select(-draw) %>%
  nest(data = I) %>%
  rename(I = data) %>%
  mutate(I = map(I, unlist)) 

# reported incidence
I_reported <- rbind(
  reported_sa %>% mutate(group = "general population"),
  reported_ch %>% mutate(group = "general population"),
  reported_ch_young %>% mutate(group = "young population")
) %>%
  left_join(nclass) %>%
  mutate(I = inci / 1e5 * n) %>%
  dplyr::select(-week, -n, -inci) %>%
  nest(data = I) %>%
  rename(I = data) %>%
  mutate(I = map(I, unlist))

# outbreak
I_outbreak <- data.frame(country = rep(names(n.class), each = n.sample),
                         n = rep(n.class, each = n.sample)) %>%
  group_by(country) %>%
  mutate(I = 1 / nclass$n[nclass$country=="South Africa"] * n) %>%
  ungroup() %>%
  dplyr::select(-n) %>%
  nest(data = I) %>%
  rename(I = data) %>%
  mutate(I = map(I, unlist))

# combine
I <- rbind(I_IFR %>% mutate(group = "IFR-based approach"),
           I_reported,
           I_outbreak %>% mutate(group = "outbreak"))

saveRDS(I, file = "data-clean/estimated-incidence-sars-cov-2.rds")


#### Summaries I ####

I %>%
  unnest() %>%
  group_by(country, group) %>%
  summarize_at("I", list(median = median, q2.5 = q2.5, q97.5 = q97.5)) %>%
  arrange(country) %>%
  mutate_if(is.numeric, round_k, 3)
  

