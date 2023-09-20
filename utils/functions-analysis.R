#### Libraries ####

library(tidyverse)
library(truncnorm)

n.sample = 3000 


#### estimating parameters ####

### rebreathed air fraction f-bar ----

f <- function(df) {
  df1 <- df %>% 
    mutate(f = ((co2 - C_o) / C_a) / 1000000,
           f.sens = ((co2 - C_o.sens) / C_a) / 1000000) %>%
    summarise(scenario = "main",
              mean = round(mean(f, na.rm = TRUE),3),
              sd = round(sd(f, na.rm = TRUE),3),
              median = round(median(f, na.rm = TRUE),3),
              quantil25 = round(quantile(f, 0.25),3),
              quantil75 = round(quantile(f, 0.75),3))
  
  df2 <- df %>%
    mutate(f = ((co2 - C_o) / C_a) / 1000000,
           f.sens = ((co2 - C_o.sens) / C_a) / 1000000) %>%
    summarise(scenario = "sensitivity",
              mean = round(mean(f.sens, na.rm = TRUE),3),
              sd = round(sd(f.sens, na.rm = TRUE),3),
              median = round(median(f.sens, na.rm = TRUE),3),
              quantil25 = round(quantile(f.sens, 0.25),3),
              quantil75 = round(quantile(f.sens, 0.75),3))
              
  bind_rows(df1, df2)
}

### I (SARS) ----


simulate_excess <- function(country, period) {
  
  # Filtering both datasets for the given country
  excess_filtered <- excess %>%
    filter(country == !!country,
           period == !!period)
  
  ifr_filtered <- ifr %>%
    filter(country == !!country,
           period == !!period)
  
  # Initial empty data frame to store results
  results <- data.frame(
    excess_draw = numeric(0),
    ifr_draq = numeric(0),
    ratio = numeric(0)
  )
  
  # Function to safely draw from truncated normal
  draw_from_truncnorm <- function(mean, sd) {
    if (mean + 1.96 * sd <= 0) {
      return(0)
    }
    return(rtruncnorm(1, a=0, mean=mean, sd=sd))
  }
  
  # Repeat process 3000 times
  for(i in 1:3000) {
    
    excess_mean <- excess_filtered$mean
    excess_sd <- excess_filtered$sd

    excess_draw <- draw_from_truncnorm(excess_mean, excess_sd)
    
    ifr_mean <- ifr_filtered $mean
    ifr_sd <- ifr_filtered $sd
    
    ifr_draw <- draw_from_truncnorm(ifr_mean, ifr_sd)
    
    # Compute the ratio and store
    results <- rbind(results, data.frame(
      excess_draw = excess_draw,
      ifr_draw = ifr_draw,
      ratio_draw = case_when(excess_draw/ifr_draw < 100000 ~ excess_draw/ifr_draw,
                        TRUE ~ 100000)
    ))
  }
  
  return(results)
}

#### Plotting ------------------------------------------------------------------

sens.df <- function(country_name) {
  
  # Validate input
  if(!country_name %in% c("South Africa", "Switzerland", "Tanzania")) {
    stop("Invalid country. Please choose one of 'South Africa', 'Switzerland', 'Tanzania'.")
  }
  
  df <- tibble(country = rep(c("South Africa", "Switzerland", "Tanzania"), each = n.sample),
               I = c(I.sa, I.ch, I.tz),
               n = rep(c(n.sa, n.ch, n.tz), each = n.sample),
               f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), each = n.sample),
               f.sens = rep(c(f_bar.sens.sa, f_bar.sens.ch, f_bar.sens.tz), each = n.sample),
               q = q.med,
               t = year) %>%
    mutate(P = case_when(
      country == {{ country_name }} ~ 1 - exp(-f.sens*q*I*t/n),
      TRUE ~ 1 - exp(-f*q*I*t/n)
    ),
    type = as.factor(case_when(
      country == {{ country_name }} ~ "600ppm",
      TRUE ~ "400ppm"
    )),
    sens = country_name)
  
  return(df)
}
