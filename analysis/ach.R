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
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.ch) / vol.ch),
         country = "Switzerland") %>% 
  summarise(mean = round(mean(ACH), 3),
            sd = round(sd(ACH),3), 
            median = round(median(ACH), 3),
            quantil25 = round(quantile(ACH, probs = 0.25), 3),
            quantil75 = round(quantile(ACH, probs = 0.75),3)) 

### Tanzania ----

ACH.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.tz) / vol.tz),
         country = "Tanzania") %>% 
  summarise(mean = round(mean(ACH), 3),
            sd = round(sd(ACH),3), 
            median = round(median(ACH), 3),
            quantil25 = round(quantile(ACH, probs = 0.25), 3),
            quantil75 = round(quantile(ACH, probs = 0.75),3))  

### South Africa ----

# Using the quantiles of the daily maxes from Switzerland

max.ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2))

# Create the ECDF for the swiss data

ecdf.ch <- ecdf(co2.ch$co2)

# Determine quantile of of the daily maxes from Switzerland

max.ch$quantile <- sapply(max.ch$co2.max, ecdf.ch)

# Using these quantiles to determine the daily maxes for South Africa

ACH.sa <- tibble(day = 1:35,
                 co2.max = quantile(co2.sa$co2, max.ch$quantile)) %>% 
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.sa) / vol.sa),
         country = "South Africa") %>% 
  summarise(mean= round(mean(ACH), 3),
            sd= round(sd(ACH),3), 
            median = round(median(ACH), 3),
            quantil25 = round(quantile(ACH, probs = 0.25), 3),
            quantil75 = round(quantile(ACH, probs = 0.75),3)) 

ACH <- bind_rows(ACH.sa, ACH.ch, ACH.tz) %>% 
  mutate(country = c("South Africa", "Switzerland", "Tanzania")) %>% 
  dplyr::select(country, everything())

write_csv(ACH, "results/co2/ACH.csv")

#### Calculation rebreathed air fraction (f) -----------------------------------

f.sa = f(co2.sa) 
f.ch = f(co2.ch)
f.tz = f(co2.tz)

f <- bind_rows(f.sa, f.ch, f.tz) %>% 
  mutate(country = rep(c("South Africa", "Switzerland", "Tanzania"), each = 2)) %>% 
  dplyr::select(country, everything())

write_csv(f, "results/co2/f.csv")

# Mean and SD of daily means of f #

mean.f_daily.ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(mean = mean(co2)) %>% 
  ungroup() %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

meanfdailych <- mean.f_daily.ch %>% 
  summarise(mean_dailymeans_f = mean(f),
            sd_dailymeans_f = sd(f))

mean.f_daily.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(mean = mean(co2)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000) 

meanfdailytz <- mean.f_daily.tz %>% 
  summarise(mean_dailymeans_f = mean(f),
            sd_dailymeans_f = sd(f))

#'different approach for South Africa, as the date isn't included in the dataset
#'described in the methods section

ecdf.ch <- ecdf(co2.ch$co2)

mean.f_daily.ch$quantile <- sapply(mean.f_daily.ch$mean, ecdf.ch)

mean.f_daily.sa <- tibble(day = 1:35,
                          mean = quantile(co2.sa$co2, mean.f_daily.ch$quantile)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

meanfdailysa <- mean.f_daily.sa %>% 
  summarise(mean_dailymeans_f = mean(f),
            sd_dailymeans_f = sd(f))

meanfdaily <- bind_rows(meanfdailysa, meanfdailych, meanfdailytz)

write_csv(meanfdaily, "results/co2/daily_means_f.csv")

# Mean and sd of daily maxes of Co2 #

max_tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  ungroup() %>% 
  summarise(mean = mean(co2.max),
            sd = sd(co2.max)) %>% 
  mutate(country = "Tanzania")

max_ch <- co2.ch %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  ungroup() %>% 
  summarise(mean = mean(co2.max),
            sd = sd(co2.max)) %>% 
  mutate(country = "Switzerland")

max_sa <- tibble(day = 1:35,
                 co2.max = quantile(co2.sa$co2, max.ch$quantile)) %>% 
  summarise(mean = mean(co2.max),
            sd = sd(co2.max)) %>% 
  mutate(country = "South Africa")

max <- bind_rows(max_sa, max_ch, max_tz)

write_csv(max, "results/co2/max.csv")
