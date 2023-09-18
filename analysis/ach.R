#### Libraries -----------------------------------------------------------------

library(tidyverse)

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
         ACH = ((3600 * Q * n.ch) / vol.ch)) %>% 
  summarise(mean=mean(ACH), 
            sd=sd(ACH), 
            median = median(ACH), 
            lower = quantile(ACH, probs = 0.25), 
            upper = quantile(ACH, probs = 0.75)) 

### Tanzania ----

ACH.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(co2.max = max(co2)) %>% 
  mutate(Q = ((0.13 * C_a * 1000000) / (co2.max - C_o)),
         ACH = ((3600 * Q * n.tz) / vol.tz)) %>% 
  summarise(mean=mean(ACH), 
            sd=sd(ACH), 
            median = median(ACH), 
            lower = quantile(ACH, probs = 0.25), 
            upper = quantile(ACH, probs = 0.75)) 

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
         ACH = ((3600 * Q * n.sa) / vol.sa)) %>% 
  summarise(mean=mean(ACH), 
            sd=sd(ACH), 
            median = median(ACH), 
            lower = quantile(ACH, probs = 0.25), 
            upper = quantile(ACH, probs = 0.75)) 

#### Calculation rebreathed air fraction (f) -----------------------------------

f.sa = f(co2.sa)
f.ch = f(co2.ch)
f.tz = f(co2.tz)
