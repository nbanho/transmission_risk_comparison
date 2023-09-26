#### Libraries -----------------------------------------------------------------

library(tidyverse)
source("utils/functions-analysis.R")

#### Data import ---------------------------------------------------------------

co2.sa <- readRDS("data-clean/co2-sa.rds") 
co2.ch <- readRDS("data-clean/co2-ch.rds") 
co2.tz <- readRDS("data-clean/co2-tz.rds") 

#### Parameters ----------------------------------------------------------------

### C_o (outdoor CO2-level) ----

C_o = 400
C_o.sens = 600

### C_a (volume fraction of CO2 added to exhaled breath during breathing) ----

C_a = ((0.0042)*60)/8

### n (number of students) ----

n.ch = 20
n.sa = 30 
n.tz = 50 

### vol (classroom volume) ----

vol.ch = 233000
vol.sa = 180000
vol.tz = 162000

#### Calculation ACH -----------------------------------------------------------

### Switzerland ----

ACH.ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  ungroup() %>%
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.ch) / vol.ch)) %>% 
  dplyr::select(-date) %>%
  summarise_all(list(mean = mean, 
                     sd = sd, 
                     med = median, 
                     q25 = function(x) quantile(x, .25), 
                     q25 = function(x) quantile(x, .75))) 

### Tanzania ----

ACH.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  ungroup() %>%
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.tz) / vol.tz)) %>%
  dplyr::select(-date) %>%
  summarise_all(list(mean = mean, 
                     sd = sd, 
                     med = median, 
                     q25 = function(x) quantile(x, .25), 
                     q25 = function(x) quantile(x, .75))) 

### South Africa ----

# Using the quantiles of the daily maxes from CH and TZ

max.ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>%
  ungroup()

max.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>%
  ungroup()

# Create the ECDF for the CH and TZ data

ecdf.ch <- ecdf(co2.ch$co2)
ecdf.tz <- ecdf(co2.tz$co2)

# Determine quantile of of the daily maxes from Switzerland

max.ch$quantile <- sapply(max.ch$co2.max, ecdf.ch)
max.tz$quantile <- sapply(max.tz$co2.max, ecdf.tz)

# Using these quantiles to determine the daily maxes for South Africa

max.sa <- quantile(co2.sa$co2, c(max.ch$quantile, max.tz$quantile))

ACH.sa <- data.frame(co2.max = c(max.sa)) %>%
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.sa) / vol.sa))  %>%
  summarise_all(list(mean = mean, 
                     sd = sd, 
                     med = median, 
                     q25 = function(x) quantile(x, .25), 
                     q25 = function(x) quantile(x, .75))) 

ACH <- bind_rows(ACH.sa, ACH.ch, ACH.tz) %>% 
  mutate(country = c("South Africa", "Switzerland", "Tanzania")) %>% 
  dplyr::select(country, everything())

ACH %>% 
  mutate_at(vars(matches("Q|ACH")), round, 2) %>%
  mutate_at(vars(matches("co2")), round, 0) %>%
  dplyr::select(country, matches("co2"), matches("Q"), matches("ACH"))

write_csv(ACH, "results/co2/ACH.csv")

#### Calculation rebreathed air fraction (f) -----------------------------------

f.sa = f(co2.sa)$mean[1]
f.ch = f(co2.ch)$mean[1]
f.tz = f(co2.tz)$mean[1]

# Mean and SD of daily means of f #

mean.f_daily.ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(C.mean = mean(co2)) %>% 
  ungroup() %>% 
  mutate(f = ((C.mean - C_o) / C_a) / 1000000)

f_daily_sum.ch <- mean.f_daily.ch %>% 
  dplyr::select(C.mean, f) %>%
  summarise_all(list(mean = mean, sd = sd))

mean.f_daily.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(C.mean = mean(co2)) %>% 
  mutate(f = ((C.mean - C_o) / C_a) / 1000000) 

f_daily_sum.tz <- mean.f_daily.tz %>% 
  dplyr::select(C.mean, f) %>%
  summarise_all(list(mean = mean, sd = sd))

#'different approach for South Africa, as the date isn't included in the dataset
#'described in the methods section

ecdf.ch <- ecdf(co2.ch$co2)
ecdf.tz <- ecdf(co2.tz$co2)

mean.f_daily.ch$quantile <- sapply(mean.f_daily.ch$C.mean, ecdf.ch)
mean.f_daily.tz$quantile <- sapply(mean.f_daily.tz$C.mean, ecdf.tz)

mean.f_daily.sa <- tibble(C.mean = quantile(co2.sa$co2, c(mean.f_daily.ch$quantile, mean.f_daily.tz$quantile))) %>% 
  mutate(f = ((C.mean - C_o) / C_a) / 1000000)

f_daily_sum.sa <- mean.f_daily.sa %>% 
  dplyr::select(C.mean, f) %>%
  summarise_all(list(mean = mean, sd = sd))

f_daily_sum <- bind_rows(f_daily_sum.sa, f_daily_sum.ch, f_daily_sum.tz) %>%
  cbind(C.mean = c(mean(co2.sa$co2), mean(co2.ch$co2), mean(co2.tz$co2))) %>%
  cbind(f = c(f.sa, f.ch, f.tz)) %>%
  mutate(country = c("South Africa", "Switzerland", "Tanzania")) %>%
  dplyr::select(country, C.mean, matches("C.mean_"), f, matches("f_"))

f_daily_sum %>% 
  mutate_at(vars(matches("f")), function(x) round(x * 100, 1)) %>%
  mutate_at(vars(matches("C.mean")), round)

write_csv(f_daily_sum, "results/co2/daily_means_f.csv")
