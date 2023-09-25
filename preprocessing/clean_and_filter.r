#### Libraries

library(haven)
library(tidyverse)

#### Data import + merging data structure ####

### SA and TZ ###

satz <- read_dta("data-raw/co2-sa-tz.dta") %>%
  filter(location == "Your school",
         grepl("class", tolower(other_location)) | other_location == "") %>%
  mutate(country = ifelse(country == 0, "South Africa", "Tanzania"),
         datetime_obj = as.POSIXct(date, format = "%Y/%m/%d %H:%M:%S"),
         date = as.Date(datetime_obj),
         time = format(datetime_obj, "%H:%M:%S"),
         time_h = hour(datetime_obj)) %>%
  dplyr::select(country, date, time, time_h, carbondioxy) %>%
  set_names(c("country", "date", "time", "time_h", "co2"))


sa <- satz %>% 
  filter(country == "South Africa",
         co2 >= 400)

tz <- satz %>% 
  filter(country == "Tanzania",
         co2 >= 400)#assume these are measurement errors

saveRDS(sa, "data-clean/co2-sa.rds")
saveRDS(tz, "data-clean/co2-tz.rds")

### CH ###

ch <- read_csv("data-raw/co2-ch.csv") %>% 
  rename(co2 = CO2) %>%
  mutate(country = "Switzerland",
         time_h = hour(time)) %>%
  filter(time_h >= 8 & time_h <= 17, #filter hours when students ar present
         co2 >= 400) %>%  #assume these are measurement errors
  dplyr::select(country, class, date, time, time_h, co2)

saveRDS(ch, "data-clean/co2-ch.rds")

#### Excess deaths file ####

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
    country = case_when(iso3c == "CHE" ~ "CH",
                        iso3c == "ZAF" ~ "SA",
                        iso3c == "TZA" ~ "TZ")) %>% 
  dplyr::select(country, date, period, estimated_daily_excess_deaths_per_100k, 
                estimated_daily_excess_deaths_ci_95_top_per_100k, 
                estimated_daily_excess_deaths_ci_95_bot_per_100k) %>% 
  mutate(weekly.per100k = estimated_daily_excess_deaths_per_100k * 7,
         weekly.per100k.low = estimated_daily_excess_deaths_ci_95_bot_per_100k * 7,
         weekly.per100k.high = estimated_daily_excess_deaths_ci_95_top_per_100k * 7) %>% 
  mutate(mean = weekly.per100k,
         sd = (weekly.per100k.high - weekly.per100k.low) / 2* 1.96) %>% 
  group_by(country, period) %>% 
  summarise(mean = mean(weekly.per100k),
            sd = sd(weekly.per100k),
            lowCI = min(weekly.per100k.low),
            highCI = min(weekly.per100k.high))

write_csv(excess, "results/prevalence/excess.csv")
