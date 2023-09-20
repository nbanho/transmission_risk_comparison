#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(truncnorm)
library(LaplacesDemon)

#### Utils ---------------------------------------------------------------------

source("utils/functions-analysis.R")
source("utils/quanta.R")
n.sample = 3000 
set.seed(1793)
country.order <- c("SA","CH", "TZ")

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

#mean girls (11-16 and 16-21), boys (11-16 and 16-21) for co2 production rate per minute.
#the breathing rate [8] comes from Rudnick paper

### f_bar (rebreathed air fraction) ----

f.sa <- f(co2.sa)
f.ch <- f(co2.ch)
f.tz <- f(co2.tz)

f_bar.sa <- f.sa$mean
f_bar.sens.sa <- f.sa$mean.sens

f_bar.ch <- f.ch$mean
f_bar.sens.ch <- f.ch$mean.sens

f_bar.tz <- f.tz$mean
f_bar.sens.tz <- f.tz$mean.sens

### n (number of students) ----

n.ch = 20
n.sa = 30 
n.tz = 50 

### t (time) ----

day = 7
week = 7*5
month = week*4
year = month*10

### q (Quanta) ----

# Distribution # 

q <- ERq(pathogen = "SARS-CoV-2", n = n.sample, pa = p_activ)

# Fixed estimates #

q.sample <- ERq(pathogen = "SARS-CoV-2", n = 1e6, pa = p_activ)
q.ann <- quantile(q.sample, c(.33, .5, .66))
scenario <- c("Low", "Median", "High")
q.med <- median(q.sample)

### I (infectious students) ----

# IFR (Lancet) #

ifr <- tibble(country = c(rep("SA", 4), rep("CH", 4), rep("TZ", 4)),
              period = rep(c("April", "July", "October", "January"), 3),
              IFR = c(0.511, 0.400, 0.348, 0.331, 1.687, 1.322, 1.164, 1.113, 0.169, 0.133, 0.116, 0.111),
              IFR.low = c(0.199, 0.173, 0.145, 0.137, 1.202, 0.943, 0.835, 0.766, 0.071, 0.060, 0.057, 0.052),
              IFR.high = c(1.275, 0.987, 0.774, 0.724, 2.782, 2.148, 1.743, 1.621, 0.339, 0.256, 0.218, 0.193)) %>% 
  mutate_at(vars(starts_with("IFR")), ~ . / 100) %>% 
  mutate(mean = IFR,
         sd = (IFR.high - IFR.low) / 2* 1.96)

# excess mortality #

excess <- read.csv("data-clean/export_country_per_100k_cumulative.csv") %>% 
  filter(iso3c %in% c("ZAF", "CHE", "TZA")) %>% 
  mutate(date = ymd(date)) %>%        # Convert date column to date format
  filter(date <= ymd("2020-12-31")) %>% 
  mutate(period = case_when(
    date >= ymd("2020-03-01") & date <= ymd("2020-04-15") ~ "April",
    date >= ymd("2020-04-16") & date <= ymd("2020-07-15") ~ "July",
    date >= ymd("2020-07-16") & date <= ymd("2020-10-15") ~ "October",
    date >= ymd("2020-10-16") & date <= ymd("2020-12-31") ~ "January",
    TRUE ~ NA_character_),
  country = case_when(iso3c == "CHE" ~ "CH",
                      iso3c == "ZAF" ~ "SA",
                      iso3c == "TZA" ~ "TZ")) %>% 
  select(country, date, period, cumulative_estimated_daily_excess_deaths_per_100k, 
         cumulative_estimated_daily_excess_deaths_ci_95_top_per_100k, 
         cumulative_estimated_daily_excess_deaths_ci_95_bot_per_100k) %>% 
  rename(weekly.per100k = cumulative_estimated_daily_excess_deaths_per_100k,
         weekly.per100k.low = cumulative_estimated_daily_excess_deaths_ci_95_bot_per_100k,
         weekly.per100k.high = cumulative_estimated_daily_excess_deaths_ci_95_top_per_100k) %>% 
  mutate(mean = weekly.per100k,
         sd = (weekly.per100k.high - weekly.per100k.low) / 2* 1.96) %>% 
  group_by(country, period) %>% 
  summarise(mean = mean(weekly.per100k),
            sd = sd(weekly.per100k))

# prevalence #

prev.april.sa <- simulate_excess("SA", "April")
prev.july.sa <- simulate_excess("SA", "July")
prev.october.sa <- simulate_excess("SA", "October")
prev.january.sa <- simulate_excess("SA", "January")

prev.april.ch <- simulate_excess("CH", "April")
prev.july.ch <- simulate_excess("CH", "July")
prev.october.ch <- simulate_excess("CH", "October")
prev.january.ch <- simulate_excess("CH", "January")

prev.april.tz <- simulate_excess("TZ", "April")
prev.july.tz <- simulate_excess("TZ", "July")
prev.october.tz <- simulate_excess("TZ", "October")
prev.january.tz <- simulate_excess("TZ", "January")

# Infectious persons in classroom per period # 

I.april.sa <- prev.april.sa$ratio_draw / 100000 * n.sa
I.july.sa <- prev.july.sa$ratio_draw / 100000 * n.sa
I.october.sa <- prev.october.sa$ratio_draw / 100000 * n.sa
I.january.sa <- prev.january.sa$ratio_draw / 100000 * n.sa

I.april.ch <- prev.april.ch$ratio_draw / 100000 * n.ch
I.july.ch <- prev.july.ch$ratio_draw / 100000 * n.ch
I.october.ch <- prev.october.ch$ratio_draw / 100000 * n.ch
I.january.ch <- prev.january.ch$ratio_draw / 100000 * n.ch

I.april.tz <- prev.april.tz$ratio_draw / 100000 * n.tz
I.july.tz <- prev.july.tz$ratio_draw / 100000 * n.tz
I.october.tz <- prev.october.tz$ratio_draw / 100000 * n.tz
I.january.tz <- prev.january.tz$ratio_draw / 100000 * n.tz

# Additional analysis #

I.sa.add = 1
I.ch.add = 1 / n.sa * n.ch
I.tz.add = 1 / n.sa * n.tz

#### Results  --------------------------------------------------------

data <- tibble(country = c(rep("SA", n.sample), rep("CH", n.sample), rep("TZ", n.sample)),
               I.april = c(I.april.sa, I.april.ch, I.april.tz),
               I.july = c(I.july.sa, I.july.ch, I.july.tz),
               I.october = c(I.october.sa, I.october.ch, I.october.tz),
               I.january = c(I.january.sa, I.january.ch, I.january.tz),
               I.add = c(rep(I.sa.add, 3000), rep(I.ch.add, 3000), rep(I.tz.add, 3000)),
               n = c(rep(n.sa, n.sample), rep(n.ch, n.sample), rep(n.tz, n.sample)),
               f = c(rep(f_bar.sa, n.sample), rep(f_bar.ch, n.sample), rep(f_bar.tz, n.sample)),
               f.sens = c(rep(f_bar.sens.sa, n.sample), rep(f_bar.sens.ch, n.sample), rep(f_bar.sens.tz, n.sample)),
               q = rep(sample.q,3),
               q100 = rep(sample.q.trunc100,3),
               q200 = rep(sample.q.trunc200,3),
               q300 = rep(sample.q.trunc300,3),
               q1 = q1,
               q2 = q2,
               q3 = q3,
               t = month,
               t.add = day,
               P.add = 1 - exp(-(f * I.add * q * t.add) / n),
               P.add100 = 1 - exp(-(f * I.add * q100 * t.add) / n),
               P.add200 = 1 - exp(-(f * I.add * q200 * t.add) / n),
               P.add300 = 1 - exp(-(f * I.add * q300 * t.add) / n),
               P.april1 = 1 - exp(-(f * I.april * q1 * t) / n),
               P.july1 = 1 - exp(-(f * I.july * q1 * t) / n),
               P.october1 = 1 - exp(-(f * I.october * q1 * t) / n),
               P.january1 = 1 - exp(-(f * I.january * q1 * t) / n),
               P.april2 = 1 - exp(-(f * I.april * q2 * t) / n),
               P.july2 = 1 - exp(-(f * I.july * q2 * t) / n),
               P.october2 = 1 - exp(-(f * I.october * q2 * t) / n),
               P.january2 = 1 - exp(-(f * I.january * q2 * t) / n),
               P.april3 = 1 - exp(-(f * I.april * q3 * t) / n),
               P.july3 = 1 - exp(-(f * I.july * q3 * t) / n),
               P.october3 = 1 - exp(-(f * I.october * q3 * t) / n),
               P.january3 = 1 - exp(-(f * I.january * q3 * t) / n))

# Main analysis #

summary.april <- data %>%
  pivot_longer(
    cols = starts_with("P.april"),
    names_to = "period",
    values_to = "P_value"
  ) %>%
  mutate(country = factor(country, levels = country.order)) %>%
  group_by(country, period) %>%
  summarise(
    ymin = quantile(P_value, 0.025),
    lower = quantile(P_value, 0.25),
    middle = quantile(P_value, 0.5),
    upper = quantile(P_value, 0.75),
    ymax = quantile(P_value, 0.975)
  )

april <- ggplot(summary.april, aes(x = country, y = middle, color = period)) + 
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), 
               stat = "identity", width = 0.5, position = position_dodge(width = 0.6), alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  ylab("") +
  theme_bw() +
  guides(fill = guide_legend(title = "Country")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel",  # This puts the title in the center of the plot area
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "April")

summary.july <- data %>%
  pivot_longer(
    cols = starts_with("P.july"),
    names_to = "period",
    values_to = "P_value"
  ) %>%
  mutate(country = factor(country, levels = country.order)) %>%
  group_by(country, period) %>%
  summarise(
    ymin = quantile(P_value, 0.025),
    lower = quantile(P_value, 0.25),
    middle = quantile(P_value, 0.5),
    upper = quantile(P_value, 0.75),
    ymax = quantile(P_value, 0.975)
  )

july <- ggplot(summary.july, aes(x = country, y = middle, color = period)) + 
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), 
               stat = "identity", width = 0.5, position = position_dodge(width = 0.6), alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  ylab("") +
  theme_bw() +
  guides(fill = guide_legend(title = "Country")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel",  # This puts the title in the center of the plot area
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "July")

summary.october <- data %>%
  pivot_longer(
    cols = starts_with("P.october"),
    names_to = "period",
    values_to = "P_value"
  ) %>%
  mutate(country = factor(country, levels = country.order)) %>%
  group_by(country, period) %>%
  summarise(
    ymin = quantile(P_value, 0.025),
    lower = quantile(P_value, 0.25),
    middle = quantile(P_value, 0.5),
    upper = quantile(P_value, 0.75),
    ymax = quantile(P_value, 0.975)
  )

october <- ggplot(summary.october, aes(x = country, y = middle, color = period)) + 
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), 
               stat = "identity", width = 0.5, position = position_dodge(width = 0.6), alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  ylab("") +
  theme_bw() +
  guides(fill = guide_legend(title = "Country")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel",  # This puts the title in the center of the plot area
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "October")

summary.january <- data %>%
  pivot_longer(
    cols = starts_with("P.january"),
    names_to = "period",
    values_to = "P_value"
  ) %>%
  mutate(country = factor(country, levels = country.order)) %>%
  group_by(country, period) %>%
  summarise(
    ymin = quantile(P_value, 0.025),
    lower = quantile(P_value, 0.25),
    middle = quantile(P_value, 0.5),
    upper = quantile(P_value, 0.75),
    ymax = quantile(P_value, 0.975)
  )

january <- ggplot(summary.january, aes(x = country, y = middle, color = period)) + 
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), 
               stat = "identity", width = 0.5, position = position_dodge(width = 0.6), alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  ylab("") +
  theme_bw() +
  guides(fill = guide_legend(title = "Country")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel",  # This puts the title in the center of the plot area
        legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "January")

months.plot <- gridExtra::grid.arrange(april, july, october, january, ncol = 2)
ggsave(plot = months.plot, filename = "results/scenarios/sars.main.png")

# Additional analysis: fixed I/n / day / quanta distribution #

summary.add <- summarize("P.add")

additional <- ggplot(summary.add, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  theme_bw() +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "",
       y = "",
       x = "")

# Additional analysis: fixed I/n / day / quanta distribution (trunc 100) #

summary.add100 <- summarize("P.add100")

additional100 <- ggplot(summary.add100, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  theme_bw() +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1))+
  labs(title = "Quanta truncated at 100 q/h",
       y = "",
       x = "")

# Additional analysis: fixed I/n / day / quanta distribution (trunc 200) #

summary.add200 <- summarize("P.add200")

additional200 <- ggplot(summary.add200, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  theme_bw() +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Quanta truncated at 200 q/h",
       y = "",
       x = "")

# Additional analysis: fixed I/n / day / quanta distribution (trunc 300) #

summary.add300 <- summarize("P.add300")

additional300 <- ggplot(summary.add300, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  theme_bw() +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1))+
  labs(title = "Quanta truncated at 300 q/h",
       y = "",
       x = "")

additional300

plot <- gridExtra::grid.arrange(additional, additional100, additional200, additional300, ncol = 2)

ggsave(plot =  plot , filename = "results/scenarios/sars.outbreak.png")

