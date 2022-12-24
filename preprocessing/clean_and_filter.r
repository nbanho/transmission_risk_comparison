library(haven)
library(tidyverse)


#### SA and TZ ####

satz <- read_dta("data-raw/co2-sa-tz.dta") %>%
  filter(location == "Your school",
         grepl("class", tolower(other_location)) | other_location == "") %>%
  mutate(country = ifelse(country == 0, "South Africa", "Tanzania"),
         date = as.character(date),
         date = as.POSIXct(date, format = "%Y/%m/%d %H:%M:%S")) %>%
  select(country, date, carbondioxy) %>%
  set_names(c("country", "date", "co2"))

saveRDS(satz, "data-clean/co2-sa-tz.rds")


#### CH ####

ch <- read_csv("data-raw/co2-ch.csv") %>%
  filter(!no_school,
         !no_class,
         intervention == "No intervention") %>%
  mutate(country = "Switzerland") %>%
  mutate(n = ifelse(class == "Study (B)", 14,
                    ifelse(class == "Study (A)", 24,
                           ifelse(school == "School 1", 14, 
                                  ifelse(school == "School 2" & class == "Study", 20, 18))))) %>%
  rename(co2 = co2ppm) %>%
  select(country, school, class, date, time, co2, n)


saveRDS(ch, "data-clean/co2-ch.rds")
