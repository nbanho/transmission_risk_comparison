#### Libraries ####

library(tidyverse)
library(truncnorm)

n.sample = 3000 


#### estimating parameters ####

### rebreathed air fraction f-bar ----

f <- function(df) {
  df %>%
    mutate(f = ((co2 - C_o) / C_a) / 1000000,
           f.sens = ((co2 - C_o.sens) / C_a) / 1000000) %>%
    summarise(mean = round(mean(f, na.rm = TRUE),3),
              sd = round(sd(f, na.rm = TRUE),3),
              median = round(median(f, na.rm = TRUE),3),
              lower = round(quantile(f, 0.25),3),
              upper = round(quantile(f, 0.75),3),
              mean.sens = round(mean(f.sens, na.rm = TRUE),3),
              sd.sens = round(sd(f.sens, na.rm = TRUE),3),
              median.sens = round(median(f.sens, na.rm = TRUE),3),
              lower.sens = round(quantile(f.sens, 0.25),3),
              upper.sens = round(quantile(f.sens, 0.75),3),)
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
