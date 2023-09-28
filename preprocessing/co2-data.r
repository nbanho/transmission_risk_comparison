#### Libraries ####

library(haven)
library(tidyverse)
source("utils/settings.R")

#### Switzerland ####

# data
ch <- read_csv("data-raw/co2-ch.csv") %>% 
  rename(co2 = CO2) %>%
  mutate(co2 = ifelse(co2 < C_o, C_o, co2)) %>%
  dplyr::select(date, class, co2)

# overall mean
ch_C.mean <- mean(ch$co2)
ch_f.mean <- ((ch_C.mean - C_o) / C_a) / 1000000
ch_f.mean.sens <- ((ch_C.mean - C_o.sens) / C_a) / 1000000

# daily mean and max
ch_C.daily <- ch %>%
  group_by(class, date) %>%
  summarize(C.mean = mean(co2),
            C.max = max(co2)) %>%
  ungroup() %>%
  mutate(f.mean = ((C.mean - C_o) / C_a) / 1000000,
         f.mean.sens = ((C.mean - C_o.sens) / C_a) / 1000000) 

ch_C.daily.sum <- ch_C.daily %>%
  summarize(across(c(C.mean, C.max, f.mean, f.mean.sens), list(mean = mean, sd = sd)))

ch_pp <- tibble(
  country = "Switzerland",
  C = list(ch$co2),
  C.mean = ch_C.mean,
  f.mean = ch_f.mean,
  f.mean.sens = ch_f.mean.sens,
  C.daily = list(ch_C.daily %>% dplyr::select(-class, -date)),
  C.daily.mean_mean = ch_C.daily.sum$C.mean_mean,
  C.daily.mean_sd = ch_C.daily.sum$C.mean_sd,
  f.daily.mean_mean = ch_C.daily.sum$f.mean_mean,
  f.daily.mean_sd = ch_C.daily.sum$f.mean_sd,
  C.daily.max_mean = ch_C.daily.sum$C.max_mean,
  C.daily.max_sd = ch_C.daily.sum$C.max_sd
)

### South Africa and Tanzania ###

# data
satz <- read_dta("data-raw/co2-sa-tz.dta") %>%
  rename(co2 = carbondioxy) %>%
  filter(location == "Your school",
         grepl("class", tolower(other_location)) | other_location == "") %>%
  mutate(country = ifelse(country == 0, "South Africa", "Tanzania"),
         co2 = ifelse(co2 < C_o, C_o, co2)) %>%
  dplyr::select(country, co2) 

# overall mean
satz.C.overall <- satz %>%
  group_by(country) %>%
  summarize(C.mean = mean(co2)) %>%
  ungroup() %>%
  mutate(f.mean = ((C.mean - C_o) / C_a) / 1000000,
         f.mean.sens = ((C.mean - C_o.sens) / C_a) / 1000000)

# daily mean and max
#' date and/or class is missing
#' --> use quantile-based approach with quantiles from empirical CDF
#' of Switzerland

ch.ecdf <- ecdf(ch_pp$C[[1]])
ch.C.daily.mean_q <- ch.ecdf(ch_pp$C.daily[[1]]$C.mean)
ch.C.daily.max_q <- ch.ecdf(ch_pp$C.daily[[1]]$C.max)

satz.C.daily <- rbind(
  data.frame(
    country = "South Africa",
    C.mean = quantile(satz$co2[satz$country=="South Africa"], ch.C.daily.mean_q),
    C.max = quantile(satz$co2[satz$country=="South Africa"], ch.C.daily.max_q)
  ),
  data.frame(
    country = "Tanzania",
    C.mean = quantile(satz$co2[satz$country=="Tanzania"], ch.C.daily.mean_q),
    C.max = quantile(satz$co2[satz$country=="Tanzania"], ch.C.daily.max_q)
  )
) %>%
  mutate(f.mean = ((C.mean - C_o) / C_a) / 1000000,
         f.mean.sens = ((C.mean - C_o.sens) / C_a) / 1000000) 

satz.C.daily.sum <- satz.C.daily %>%
  group_by(country) %>%
  summarize(across(c(C.mean, C.max, f.mean, f.mean.sens), list(mean = mean, sd = sd)))

# South Africa
sa_pp <- tibble(
  country = "South Africa",
  C = list(satz$co2[satz$country=="South Africa"]),
  C.mean = satz.C.overall$C.mean[satz.C.overall$country=="South Africa"],
  f.mean = satz.C.overall$f.mean[satz.C.overall$country=="South Africa"],
  f.mean.sens = satz.C.overall$f.mean.sens[satz.C.overall$country=="South Africa"],
  C.daily = list(satz.C.daily %>% filter(country == "South Africa") %>% dplyr::select(-country)),
  C.daily.mean_mean = satz.C.daily.sum$C.mean_mean[satz.C.daily.sum$country=="South Africa"],
  C.daily.mean_sd = satz.C.daily.sum$C.mean_sd[satz.C.daily.sum$country=="South Africa"],
  f.daily.mean_mean = satz.C.daily.sum$f.mean_mean[satz.C.daily.sum$country=="South Africa"],
  f.daily.mean_sd = satz.C.daily.sum$f.mean_sd[satz.C.daily.sum$country=="South Africa"],
  C.daily.max_mean = satz.C.daily.sum$C.max_mean[satz.C.daily.sum$country=="South Africa"],
  C.daily.max_sd = satz.C.daily.sum$C.max_sd[satz.C.daily.sum$country=="South Africa"]
)

# Tanzania
tz_pp <- tibble(
  country = "Tanzania",
  C = list(satz$co2[satz$country=="Tanzania"]),
  C.mean = satz.C.overall$C.mean[satz.C.overall$country=="Tanzania"],
  f.mean = satz.C.overall$f.mean[satz.C.overall$country=="Tanzania"],
  f.mean.sens = satz.C.overall$f.mean.sens[satz.C.overall$country=="Tanzania"],
  C.daily = list(satz.C.daily %>% filter(country == "Tanzania") %>% dplyr::select(-country)),
  C.daily.mean_mean = satz.C.daily.sum$C.mean_mean[satz.C.daily.sum$country=="Tanzania"],
  C.daily.mean_sd = satz.C.daily.sum$C.mean_sd[satz.C.daily.sum$country=="Tanzania"],
  f.daily.mean_mean = satz.C.daily.sum$f.mean_mean[satz.C.daily.sum$country=="Tanzania"],
  f.daily.mean_sd = satz.C.daily.sum$f.mean_sd[satz.C.daily.sum$country=="Tanzania"],
  C.daily.max_mean = satz.C.daily.sum$C.max_mean[satz.C.daily.sum$country=="Tanzania"],
  C.daily.max_sd = satz.C.daily.sum$C.max_sd[satz.C.daily.sum$country=="Tanzania"]
)

#### Save ####

co2 <- rbind(ch_pp, sa_pp, tz_pp) %>%
  mutate(country = factor(country, levels = ctry_ord)) %>%
  arrange(country)

saveRDS(co2, "data-clean/co2.rds")


#### Reporting ####

co2 %>% 
  dplyr::select(-C, -C.daily) %>%
  mutate_at(vars(matches("C[[:punct:]]")), round_k, 0) %>%
  mutate_at(vars(matches("f")), function(x) round_k(x * 100, 1)) 

