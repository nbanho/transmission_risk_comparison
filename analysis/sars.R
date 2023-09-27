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

#' mean girls (11-16 and 16-21), boys (11-16 and 16-21) for co2 production rate per minute.
#' the breathing rate [8] comes from Rudnick paper

### f_bar (rebreathed air fraction) ----

f.sa <- f(co2.sa)
f.ch <- f(co2.ch)
f.tz <- f(co2.tz)

# Main analysis #

f_bar.sa <- f.sa %>% 
  filter(scenario == "main") %>% 
  pull(mean)

f_bar.ch <- f.ch %>% 
  filter(scenario == "main") %>% 
  pull(mean)

f_bar.tz <- f.tz %>% 
  filter(scenario == "main") %>% 
  pull(mean)

# Additional analysis (outbreak) #

mean.f_daily.ch <- co2.ch %>% 
  group_by(date, class) %>% 
  summarise(mean = mean(co2)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

#'different approach for South Africa and Tanzania, as the date isn't included in the dataset
#'described in the methods section

ecdf.ch <- ecdf(co2.ch$co2)

mean.f_daily.ch$quantile <- sapply(mean.f_daily.ch$mean, ecdf.ch)

mean.f_daily.tz <- tibble(C.mean = quantile(co2.tz$co2, mean.f_daily.ch$quantile)) %>% 
  mutate(f = ((C.mean - C_o) / C_a) / 1000000)

mean.f_daily.sa <- tibble(C.mean = quantile(co2.sa$co2, mean.f_daily.ch$quantile)) %>% 
  mutate(f = ((C.mean - C_o) / C_a) / 1000000)

f.sample.sa <- sample(mean.f_daily.sa$f , size = 9000, replace = TRUE)
f.sample.ch <- sample(mean.f_daily.ch$f , size = 9000, replace = TRUE)
f.sample.tz <- sample(mean.f_daily.tz$f , size = 9000, replace = TRUE)

# Sensitivity analysis (outdoor CO2) #

f_bar.sens.sa <- f.sa %>% 
  filter(scenario == "sensitivity") %>% 
  pull(mean)

f_bar.sens.ch <- f.ch %>% 
  filter(scenario == "sensitivity") %>% 
  pull(mean)

f_bar.sens.tz <- f.tz %>% 
  filter(scenario == "sensitivity") %>% 
  pull(mean)

### n (number of students) ----

n.ch = 20
n.sa = 30 
n.tz = 50 

### t (time) ----

day = 919 / 365
week = 919 / 52
month = 919 / 12
year = 919

### q (Quanta) ----

# Distribution # 

q.sars.low <- ERq("SARS-CoV-2", n.sample, pa = p_low)
q.sars.med <- ERq("SARS-CoV-2", n.sample, pa = p_med)
q.sars.high <- ERq("SARS-CoV-2", n.sample, pa = p_high)

q_long <- tibble(
  country = rep(c("South Africa", "Switzerland", "Tanzania"), each = n.sample * 3),
  scenario = rep(rep(c("Low", "Medium", "High"), each = n.sample), times = 3),
  q = rep(c(unlist(q.sars.low), unlist(q.sars.med), unlist(q.sars.high)),3)
)

### I (infectious students) ----

# prevalence #

prev.sa <- simulate_excess("SA")
prev.ch <- simulate_excess("CH")
prev.tz <- simulate_excess("TZ")

# calculating I #

I.sa <- prev.sa %>%
  group_by(draw) %>% 
  summarise(mean.prev = mean(draw_prev_per100k)) %>% 
  mutate(I_weekly = mean.prev / 100000 * n.sa,
         country = "SA")

I.ch <- prev.ch %>%
  group_by(draw) %>% 
  summarise(mean.prev = mean(draw_prev_per100k)) %>% 
  mutate(I_weekly = mean.prev / 100000 * n.ch,
         country = "CH")

I.tz <- prev.tz %>%
  group_by(draw) %>% 
  summarise(mean.prev = mean(draw_prev_per100k)) %>% 
  mutate(I_weekly = mean.prev / 100000 * n.tz,
         country = "TZ")

I.results <- bind_rows(I.sa, I.ch, I.tz) %>% 
  group_by(country) %>% 
  summarise(lowCI_weekly = quantile(I_weekly, 0.025),
            median_weekly = median(I_weekly),
            highCI_weekly = quantile(I_weekly, 0.975))

I.results %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(country = factor(country, levels = c("SA", "CH", "TZ"))) %>%
  arrange(country) %>%
  dplyr::select(country, median_weekly, lowCI_weekly, highCI_weekly)

# Additional analysis #

I.sa.add = 1
I.ch.add = 1 / n.sa * n.ch
I.tz.add = 1 / n.sa * n.tz

# Sensitivity analysis (reported cases) #

reported <- readRDS("data-clean/reported_weekly.rds")
reported_ch <- readRDS("data-clean/reported_weekly_ch.rds")
reported_young <- readRDS("data-clean/reported_weekly_young.rds")                        

I.sens.sa <- tibble(I = c(sample(reported %>% 
                              filter(Country == "South Africa") %>% 
                              group_by(week_start) %>% 
                              summarise(cases_weekly = sum(New_cases)) %>% 
                              pull(cases_weekly), size = 3000, replace = TRUE) / 59390000 * n.sa,
                           I.sa$I_weekly,
                          rep(0, 3000)), scenario = rep(c("reported", "IFR", "reported young"), each = 3000))

I.sens.ch <- tibble(I = c(sample(reported_ch %>%
                                   pull(weekly_cases), size = 3000, replace = TRUE)/ 8700000 * n.ch,
                          I.ch$I_weekly,
                          sample(reported_young$weekly_cases, size = 3000, replace = TRUE) / 1755000 * n.ch), 
                    scenario = rep(c("reported", "IFR", "reported young"), each = 3000))

I.sens.tz <- tibble(I = c(I.tz$I_weekly, rep(0, 6000)),
                    scenario = rep(c("IFR", "reported", "reported young"), each = 3000))

I.sens_sum <- rbind(
  I.sens.sa %>% dplyr::filter(scenario == "reported") %>% mutate(country = "South Africa"),
  I.sens.ch %>% dplyr::filter(scenario != "IFR") %>% mutate(country = "Switzerland")
) %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania"))) %>%
  group_by(country, scenario) %>%
  summarize(medianI = median(I),
            Q2.5 = quantile(I, 0.025),
            Q97.5 = quantile(I, .975))

I.sens_sum %>% mutate_if(is.numeric, round, 2)

#### Results  ------------------------------------------------------------------

# Main analysis #

df.main <- tibble(
  country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                each = n.sample * 3),
  I = rep(c(I.sa$I_weekly, I.ch$I_weekly, I.tz$I_weekly), each = 3),
  n = rep(c(n.sa, n.ch, n.tz), 
          each = n.sample * 3),
  f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
          each = n.sample * 3),
  t = month,
  q = q_long$q,
  sc = q_long$scenario) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))

plottingMtb(df.main) +
  ylab("Monthly risk of infection (%, square-root scale)")

ggsave("results/SARS-CoV-2/sars-main.png", width = 16 / cm(1), height = 10 / cm(1))  

results.main <- results(df.main)

results.main

write_csv(results.main, "results/SARS-CoV-2/main.csv")

# Additional analysis: fixed I/n / day / sampling daily means of f #

df.add <- tibble(
  country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                each = n.sample * 3),
  I = rep(c(I.sa.add, I.ch.add, I.tz.add),
          each = n.sample * 3),
  n = rep(c(n.sa, n.ch, n.tz), 
          each = n.sample * 3),
  f = c(f.sample.sa, f.sample.ch, f.sample.tz),
  t = week,
  q = q_long$q,
  sc = q_long$scenario) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n)) 

plottingMtb(df.add) +
  labs(y = "Weekly risk of infection (%, square-root scale)")

ggsave("results/SARS-CoV-2/sars.add.png", width = 16 / cm(1), height = 10 / cm(1))  

results.add <- results(df.add)

results.add

write_csv(results.add, "results/SARS-CoV-2/add.csv")

# Sensitivity analysis (Incidence): Comparing IFR-approach / reported incidence all age / reported incidence young people

df.sens.prev <- tibble(
  country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                each = n.sample * 3),
  I = rep(c(I.sens.sa$I, I.sens.ch$I, I.sens.tz$I)),
  sc = rep(c(I.sens.sa$scenario, I.sens.ch$scenario, I.sens.tz$scenario)),
  n = rep(c(n.sa, n.ch, n.tz), 
          each = n.sample * 3),
  f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
          each = n.sample * 3),
  t = month,
  q = rep(q.sars.med, 3*3)) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n)) %>%
  mutate(sc = factor(sc, levels = c("reported", "reported young", "IFR")))

ggplot(mapping = aes(x = sc, y = P, fill = country)) +
  geom_errorbar(data = df.sens.prev %>% 
                  dplyr::select(country, sc, P) %>% 
                  group_by(country, sc) %>% 
                  median_qi(),
                mapping = aes(ymin = .lower, ymax = .upper, color = country),
                position = position_dodge(width = .5),
                width = .3) +
  geom_boxplot(data = df.sens.prev, 
               position = position_dodge(width = .5),
               outlier.shape = NA, coef = 0, width = 0.3) +
  stat_summary(data = df.sens.prev,
               mapping = aes(color = country), 
               geom = "point", 
               fun = function(y) ifelse(median(y) == 0, NA, median(y)), 
               position = position_dodge2(width = .5), 
               size = 2, shape = 23, fill = "white")+
  scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
               limits = c(-1,1), breaks = seq(0, 1, .2)^2) +
  scale_x_discrete(labels = c("Reported incidence\n(general population)", "Reported incidence\n(age group 10-20y)", "IFR-based approach"), expand = expansion(add = c(.33, .33))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(y = "Monthly risk of infection (%, square-root scale)", 
       x = "",
       color = '',
       fill = '') +
  theme_custom() +
  theme(legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") 

ggsave("results/SARS-CoV-2/sars.sens_prev.png", width = 16 / cm(1), height = 10 / cm(1))  

results.sens.prev <- results(df.sens.prev)

results.sens.prev

write_csv(results.add, "results/SARS-CoV-2/sens_prev.csv")

# Sensitivity analysis (Outdoor CO2): Comparing 600 ppm vs 400 ppm #

df.sens.co2.sa <- sens.df.sars("South Africa")
df.sens.co2.ch <- sens.df.sars("Switzerland")
df.sens.co2.tz <- sens.df.sars("Tanzania")

df.sens.co2 <- bind_rows(df.sens.co2.sa, df.sens.co2.ch, df.sens.co2.tz) 

df.sens.co2 %>% 
  mutate(sc = factor(sc, levels = c("South Africa", "Switzerland", "Tanzania"))) %>% 
  ggplot(mapping = aes(x = sc, y = P, fill = country, shape = type)) +
  geom_errorbar(data = df.sens.co2 %>% 
                  dplyr::select(country, sc, P, type) %>% 
                  group_by(country, sc, type) %>% 
                  median_qi(),
                mapping = aes(ymin = .lower, ymax = .upper, color = country),
                position = position_dodge(width = .5),
                width = .3) +
  geom_boxplot(data = df.sens.co2, 
               position = position_dodge(width = .5),
               outlier.shape = NA, coef = 0, width = 0.3) +
  stat_summary(data = df.sens.co2,
               mapping = aes(x = sc, y = P, group = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 4, fill = "white", color = "black") +
  scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
               limits = c(0,1), breaks = seq(0, 1, .2)^2) +
  scale_x_discrete(labels = c(expression(C^o*"=600ppm in South Africa"), 
                              expression(C^o*"=600ppm in Switzerland"), 
                              expression(C^o*"=600ppm in Tanzania")), 
                   expand = expansion(add = c(.33, .33))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  scale_shape_manual(values=c(23, 17), name = expression("Outdoor "*CO[2]*"-Level")) +  
  labs(y = "Annual risk of infection (%, square-root scale)", 
       x = "",
       color = '',
       fill = '') +
  theme_custom() +
  theme(legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") 

ggsave("results/SARS-CoV-2/mtb-sens.co2.png", width = 16 / cm(1), height = 10 / cm(1))  

results.sens.co2 <- results(df.sens.co2)

write_csv(results.sens.co2, "results/SARS-CoV-2/sens_co2.csv")
