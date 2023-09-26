#### Libraries ####

library(tidyverse)
library(truncnorm)
library(purrr)

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

# IFR (Lancet) #

excess <- readRDS("data-clean/excess.rds")

ifr <- tibble(country = c(rep("SA", 4), rep("CH", 4), rep("TZ", 4)),
              period = rep(c("April", "July", "October", "January"), 3),
              IFR = c(0.511, 0.400, 0.348, 0.331, 1.687, 1.322, 1.164, 1.113, 0.169, 0.133, 0.116, 0.111),
              IFR.low = c(0.199, 0.173, 0.145, 0.137, 1.202, 0.943, 0.835, 0.766, 0.071, 0.060, 0.057, 0.052),
              IFR.high = c(1.275, 0.987, 0.774, 0.724, 2.782, 2.148, 1.743, 1.621, 0.339, 0.256, 0.218, 0.193)) %>% 
  mutate_at(vars(starts_with("IFR")), ~ . / 100) %>% 
  mutate(mean = IFR,
         sd = (IFR.high - IFR.low) / 2* 1.96)

# excess mortality #

simulate_excess <- function(country) {
  
  unique_periods <- unique(excess$period)
  
  # Function to draw from truncated normal
  draw_from_truncnorm <- function(mean, sd) {
    return(rtruncnorm(1, a=0, mean=mean, sd=sd))
  }
  
  # Create an empty results dataframe
  results <- data.frame(
    draw_excess = numeric(0),
    draw_ifr = numeric(0),
    draw_prev_per100k = numeric(0),
    period = character(0),
    draw = integer(0)
  )
  
  # Repeat process 3000 times
  for(i in 1:3000) {
    
    # Simulate for each period
    for (period1 in unique_periods) {
      
      # Filtering both datasets for the given country and period
      excess_filtered <- excess %>%
        filter(country == !!country,
               period == period1)
      
      ifr_filtered <- ifr %>%
        filter(country == !!country,
               period == period1)
      
      excess_mean <- excess_filtered$mean
      excess_sd <- excess_filtered$sd
      excess_draw <- draw_from_truncnorm(excess_mean, excess_sd)
      
      ifr_mean <- ifr_filtered$mean
      ifr_sd <- ifr_filtered$sd
      ifr_draw <- draw_from_truncnorm(ifr_mean, ifr_sd)
      
      # Compute the ratio and store, including the draw number
      results <- rbind(results, data.frame(
        draw_excess = excess_draw,
        draw_ifr = ifr_draw,
        draw_prev_per100k = pmin(excess_draw/ifr_draw, 100000),
        period = period1,
        draw = i
      ) %>% 
       dplyr:: select(draw, period, draw_excess, draw_ifr, draw_prev_per100k))
    }
  }
  
  return(results)
}

#### Plotting ------------------------------------------------------------------

## Mtb ##

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
               q = rep(q.mtb.med, 3),
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

## SARS-CoV-2 ##

sens.df.sars <- function(country_name) {
  
  # Validate input
  if(!country_name %in% c("South Africa", "Switzerland", "Tanzania")) {
    stop("Invalid country. Please choose one of 'South Africa', 'Switzerland', 'Tanzania'.")
  }
  
  df <- tibble(country = rep(c("South Africa", "Switzerland", "Tanzania"), each = n.sample),
               I = c(I.sa$I_weekly, I.ch$I_weekly, I.tz$I_weekly),
               n = rep(c(n.sa, n.ch, n.tz), each = n.sample),
               f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), each = n.sample),
               f.sens = rep(c(f_bar.sens.sa, f_bar.sens.ch, f_bar.sens.tz), each = n.sample),
               q = rep(q.sars.med, 3),
               t = month) %>%
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

# main and additional plot #

plottingMtb <- function(df) {
  
ggplot(mapping = aes(x = factor(sc, levels = c("Low", "Medium", "High")), y = P, fill = country)) +
  geom_errorbar(data = df %>% 
                  dplyr::select(country, sc, P) %>% 
                  group_by(country, sc) %>% 
                  median_qi(),
                mapping = aes(ymin = .lower, ymax = .upper, color = country),
                position = position_dodge(width = .5),
                width = .3) +
  geom_boxplot(data = df, 
               position = position_dodge(width = .5),
               outlier.shape = NA, coef = 0, width = 0.3) +
  stat_summary(data = df,
               mapping = aes(x = sc, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
               limits = c(0,1), breaks = seq(0, 1, .2)^2) +
  scale_x_discrete(labels = c("Low", "Medium", "High"), expand = expansion(add = c(.33, .33))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(y = "Annual risk of infection (%, square-root scale)", 
       x = "Activity level in the classroom",
       color = '',
       fill = '') +
  theme_custom() +
  theme(legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") 
  
}

plottingSars <- function(df) {
  
  ggplot(mapping = aes(x = factor(sc, levels = c("Low", "Medium", "High")), y = P, fill = country)) +
    geom_errorbar(data = df %>% 
                    dplyr::select(country, sc, P) %>% 
                    group_by(country, sc) %>% 
                    median_qi(),
                  mapping = aes(ymin = .lower, ymax = .upper, color = country),
                  position = position_dodge(width = .5),
                  width = .3) +
    geom_boxplot(data = df, 
                 position = position_dodge(width = .5),
                 outlier.shape = NA, coef = 0, width = 0.3) +
    stat_summary(data = df,
                 mapping = aes(x = sc, y = P, color = country), geom = "point", fun = "median", 
                 position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
    scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
                 limits = c(0,1), breaks = seq(0, 1, .2)) +
    scale_x_discrete(labels = c("Low", "Medium", "High"), expand = expansion(add = c(.33, .33))) +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2")) +
    labs(y = "Monthly risk of infection (%)", 
         x = "Activity level in the classroom",
         color = '',
         fill = '') +
    theme_custom() +
    theme(legend.position = "top",
          plot.title.position = "plot",
          legend.box = "vertical") 
  
}

# results dataframe #
results <- function(df){
  
  df %>% 
    group_by(country) %>% 
    dplyr::summarise(lowerCI = round(quantile(P, 0.025),3),
                     median = round(median(P), 3),
                     upperCI = round(quantile(P, 0.975),3))
  
}
