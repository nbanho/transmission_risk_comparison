#### TO-DO ####

#outbreak scenario: sampling from daily f means instead of fixed f

#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(truncnorm)
library(LaplacesDemon)
library(tidybayes)
library(wesanderson)

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
  group_by(date) %>% 
  summarise(mean = mean(co2)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

mean.f_daily.tz <- co2.tz %>% 
  group_by(date) %>% 
  summarise(mean = mean(co2)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

#'different approach for South Africa, as the date isn't included in the dataset
#'described in the methods section

ecdf.ch <- ecdf(co2.ch$co2)

mean.f_daily.ch$quantile <- sapply(mean.f_daily.ch$mean, ecdf.ch)

mean.f_daily.sa <- tibble(day = 1:35,
                 mean = quantile(co2.sa$co2, mean.f_daily.ch$quantile)) %>% 
  mutate(f = ((mean - C_o) / C_a) / 1000000)

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

day = 7
week = 7*5
month = week*4
year = 919

### q (Quanta) ----

# Distribution # 

q.mtb.low <- ERq("Mtb", n.sample, pa = p_low)
q.mtb.med <- ERq("Mtb", n.sample, pa = p_med)
q.mtb.high <- ERq("Mtb", n.sample, pa = p_high)

q_long <- tibble(
  country = rep(c("South Africa", "Switzerland", "Tanzania"), each = n.sample * 3),
  scenario = rep(rep(c("Low", "Medium", "High"), each = n.sample), times = 3),
  q = rep(c(unlist(q.mtb.low), unlist(q.mtb.med), unlist(q.mtb.high)),3)
)

### I (infectious students) ----

# Main analysis #

prev.sa <- rtruncnorm(n.sample, mean = 432, sd = (632 - 232) / (2 * qnorm(0.975)))
prev.ch <- rtruncnorm(n.sample, a = 0, mean = 12, sd = 3.7)
prev.tz <- rtruncnorm(n.sample, a = 0, mean = 42, sd = (73 - 11) / (2 * qnorm(0.975)))

I.sa = prev.sa / 100000 * n.sa
I.ch = prev.ch / 100000 * n.ch
I.tz = prev.tz / 100000 * n.tz

# Sensitivity analysis (prevalence) #

prev.sa.sens <- rtruncnorm(n.sample, a = 0, mean = 852, sd = ((1026 - 679) / (2 * qnorm(0.975))))
prev.ch.sens <- rtruncnorm(n.sample, a = 0, mean = 56.9, sd = 6.4)
prev.tz.sens <- rtruncnorm(n.sample, a = 0, mean = 293, sd = ((358 - 228) / (2 * qnorm(0.975))))

I.sa.sens = prev.sa.sens / 100000 * n.sa
I.ch.sens = prev.ch.sens / 100000 * n.ch
I.tz.sens = prev.tz.sens / 100000 * n.tz

# Additional analysis (outbreak) #

I.sa.add = 1
I.ch.add = 1 / n.sa * n.ch
I.tz.add = 1 / n.sa * n.tz

#### Results  ------------------------------------------------------------------

# Main analysis #

df.main <- tibble(
                country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                                each = n.sample * 3),
                I = rep(c(I.sa, I.ch, I.tz), each = 3),
                n = rep(c(n.sa, n.ch, n.tz), 
                        each = n.sample * 3),
                f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
                        each = n.sample * 3),
                t = year,
                q = q_long$q,
                sc = q_long$scenario) %>% 
            mutate(P = 1 - exp(-f*q*I*t/n))
  
plottingMtb(df.main)

ggsave("results/Mtb/mtb-main.png", width = 12 / cm(1), height = 8 / cm(1))  

results.main <- results(df.main)

write_csv(results.main, "results/Mtb/main.csv")
# Additional analysis (outbreak): fixed I/n / day / sampling daily means of f #

df.add <- tibble(
                  country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                               each = n.sample * 3),
                  I = rep(c(I.sa.add, I.ch.add, I.tz.add),
                          each = n.sample* 3),
                  n = rep(c(n.sa, n.ch, n.tz), 
                          each = n.sample * 3),
                  f = c(f.sample.sa, f.sample.ch, f.sample.tz),
                  t = week,
                  q = q_long$q,
                  sc = q_long$scenario) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))

plottingMtb(df.add) +
  ylab("Weekly risk of infection (%, square root scale)")

ggsave("results/Mtb/mtb-add.png", width = 12 / cm(1), height = 8 / cm(1))  

results.add <- results(df.add)

write_csv(results.add, "results/Mtb/add.csv")

# Sensitivity analysis (prevalence): Using the non-age-specific prevalence data #

df.sens.prev <- tibble(country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                               each = n.sample),
                 I = c(I.sa.sens, I.ch.sens, I.tz.sens),
                 n = rep(c(n.sa, n.ch, n.tz), 
                         each = n.sample),
                 f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
                         each = n.sample),
                 f.sens = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
                              each = n.sample),
                 q = rep(q.mtb.med, 3),
                 t = year) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))

ggplot(mapping = aes(x = factor(country, levels = c("South Africa", "Switzerland", "Tanzania")), y = P, fill = country)) +
  geom_errorbar(data = df.sens.prev %>% 
                  dplyr::select(country, P) %>% 
                  group_by(country) %>% 
                  median_qi(),
                mapping = aes(ymin = .lower, ymax = .upper, color = country),
                position = position_dodge(width = .5),
                width = .3) +
  geom_boxplot(data = df, 
               position = position_dodge(width = .5),
               outlier.shape = NA, coef = 0, width = 0.3) +
  stat_summary(data = df,
               mapping = aes(x = country, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
               limits = c(0,1), breaks = seq(0, 1, .2)^2) +
  scale_x_discrete(labels = c("South Africa", "Switzerland", "Tanzania"), expand = expansion(add = c(.33, .33))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(y = "Annual risk of infection (%, square-root scale)", 
       x = "",
       color = '',
       fill = '') +
  theme_custom() +
  theme(legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") 

ggsave("results/Mtb/mtb-sens.prev.png", width = 12 / cm(1), height = 8 / cm(1))  

results.sens.prev <- results(df.sens.prev)

write_csv(results.sens.prev, "results/Mtb/sens.prev.csv")

# Sensitivity analysis (Outdoor CO2): Comparing 600 ppm vs 400 ppm #

df.sens.co2.sa <- sens.df("South Africa")
df.sens.co2.ch <- sens.df("Switzerland")
df.sens.co2.tz <- sens.df("Tanzania")

df.sens.co2 <- bind_rows(df.sens.co2.sa, df.sens.co2.ch, df.sens.co2.tz) 

df.sens.co2 %>% 
  mutate(sens = factor(sens, levels = c("South Africa", "Switzerland", "Tanzania"))) %>% 
  ggplot(mapping = aes(x = sens, y = P, fill = country, shape = type)) +
    geom_errorbar(data = df.sens.co2 %>% 
                    dplyr::select(country, sens, P, type) %>% 
                    group_by(country, sens, type) %>% 
                    median_qi(),
                  mapping = aes(ymin = .lower, ymax = .upper, color = country),
                  position = position_dodge(width = .5),
                  width = .3) +
    geom_boxplot(data = df.sens.co2, 
                 position = position_dodge(width = .5),
                 outlier.shape = NA, coef = 0, width = 0.3) +
    stat_summary(data = df.sens.co2,
                 mapping = aes(x = sens, y = P, color = country), geom = "point", fun = "median", 
                 position = position_dodge2(width = .5), size = 2, fill = "white") +
    scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
                 limits = c(0,1), breaks = seq(0, 1, .2)^2) +
    scale_x_discrete(labels = c("South Africa", "Switzerland", "Tanzania"), expand = expansion(add = c(.33, .33))) +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2")) +
    scale_shape_manual(values=c(23, 24), name = expression("Outdoor" ~ CO[2]- ~ "Level")) +  
    labs(y = "Annual risk of infection (%, square-root scale)", 
         x = "",
         color = '',
         fill = '') +
    theme_custom() +
    theme(legend.position = "top",
          plot.title.position = "plot",
          legend.box = "vertical") 

ggsave("results/Mtb/mtb-sens.co2.png", width = 12 / cm(1), height = 8 / cm(1))  

results.sens.co2 <- results(df.sens.co2)

write_csv(results.sens.co2, "results/Mtb/sens.co2.csv")

