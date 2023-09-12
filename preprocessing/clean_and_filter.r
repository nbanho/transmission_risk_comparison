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
