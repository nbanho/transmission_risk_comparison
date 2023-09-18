#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(truncnorm)
library(LaplacesDemon)
library(tidybayes)

#### Utils ---------------------------------------------------------------------

source("utils/functions-analysis.R")
source("utils/quanta.R")
source("utils/plotting.R")
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

f_bar.sa <- f.sa$mean
f_bar.ch <- f.ch$mean
f_bar.tz <- f.tz$mean

# Sensitivity analysis (outdoor CO2) #

f_bar.sens.sa <- f.sa$mean.sens
f_bar.sens.ch <- f.ch$mean.sens
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

q <- ERq(pathogen = "Mtb", n = n.sample, pa = p_activ)

# Fixed estimates #

q.sample <- ERq(pathogen = "Mtb", n = 1e6, pa = p_activ)
q.ann <- quantile(q.sample, c(.33, .5, .66))
scenario <- c("Low", "Median", "High")
q.med <- median(q.sample)

### I (infectious students) ----

# Main analysis #

prev.sa <- rtruncnorm(n.sample, mean = 432, sd = (632 - 232) / (2 * qnorm(0.975)))
prev.ch <- rtruncnorm(n.sample, a = 0, mean = 12, sd = 3.7)
prev.tz <- rtruncnorm(n.sample, a = 0, mean = 42, sd = (73 - 11) / (2 * qnorm(0.975)))

I.sa = prev.sa / 100000 * n.sa
I.ch = prev.ch / 100000 * n.ch
I.tz = prev.tz / 100000 * n.tz

# Sensitivity analysis (prevalence) #

prev.sa.sens <- rtruncnorm(n.sample, a = 0, mean = 852, sd = ((1026 - 679) / 2 * qnorm(0.975)))
prev.ch.sens <- rtruncnorm(n.sample, a = 0, mean = 56.9, sd = 6.4)
prev.tz.sens <- rtruncnorm(n.sample, a = 0, mean = 293, sd = ((358 - 228) / 2 * qnorm(0.975)))

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
                                each = n.sample),
                I = c(I.sa, I.ch, I.tz),
                n = rep(c(n.sa, n.ch, n.tz), 
                        each = n.sample),
                f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
                        each = n.sample),
                q = c(rep(list(q.ann), 
                          n.sample*3)),
                sc = c(rep(list(scenario), 
                           n.sample*3)),
                t = year) %>% 
  unnest(c(q, sc)) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))
  
df.main %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania")),
         sc = factor(sc, levels = c("Low", "Median", "High"))) %>%
  ggplot(aes(x = sc, y = P)) +
  stat_interval(aes(color = country, color_ramp = after_stat(rev(.width))), position = position_dodge(width = .5)) +
  stat_summary(aes(x = sc, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0.05)), limits = c(0,1)) +
  scale_x_discrete(labels = c(expression(atop("Low", q*' = 0.35'~'h'^-1)), 
                              expression(atop("Median", q*' = 1.53'~'h'^-1)),
                              expression(atop("High", q*' = 6.26'~'h'^-1)))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  ggdist::scale_color_ramp_continuous() +
  labs(y = "Risk of infection (%)", color_ramp = "Interval", color = "") +
  theme_custom() +
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") +
  guides(color_ramp = "none", 
         color = guide_legend(order = 1))

ggsave("results/mtb-main.png", width = 12 / cm(1), height = 8 / cm(1))  

# Additional analysis (outbreak): fixed I/n / day / quanta distribution #

df.add <- tibble(country = rep(c("South Africa", "Switzerland", "Tanzania"), 
                               each = n.sample),
                  I = rep(c(I.sa.add, I.ch.add, I.tz.add),
                          each = n.sample),
                  n = rep(c(n.sa, n.ch, n.tz), 
                          each = n.sample),
                  f = rep(c(f_bar.sa, f_bar.ch, f_bar.tz), 
                          each = n.sample),
                  q = rep(q, 3),
                  t = week) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))

df.add %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania"))) %>%
  ggplot(aes(x = country, y = P)) +
  stat_interval(aes(color = country, color_ramp = after_stat(rev(.width))), position = position_dodge(width = .5)) +
  stat_summary(aes(x = country, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0.05)), limits = c(0,1)) +
  scale_x_discrete(labels = c(expression(atop("South Africa")), 
                              expression(atop("Switzerland")),
                              expression(atop("Tanzania")))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  ggdist::scale_color_ramp_continuous() +
  labs(y = "Risk of infection (%)", color_ramp = "Interval", color = "") +
  theme_custom() +
  theme(axis.title.x = element_blank()) +
  guides(color_ramp = "none", 
         color = "none")

ggsave("results/mtb-add.png", width = 12 / cm(1), height = 8 / cm(1))  

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
                 q = q.med,
                 t = year) %>% 
  mutate(P = 1 - exp(-f*q*I*t/n))

df.sens.prev %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania"))) %>%
  ggplot(aes(x = country, y = P)) +
  stat_interval(aes(color = country, color_ramp = after_stat(rev(.width))), position = position_dodge(width = .5)) +
  stat_summary(aes(x = country, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0.05)), limits = c(0,1)) +
  scale_x_discrete(labels = c(expression(atop("South Africa")), 
                              expression(atop("Switzerland")),
                              expression(atop("Tanzania")))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  ggdist::scale_color_ramp_continuous() +
  labs(y = "Risk of infection (%)", color_ramp = "Interval", color = "") +
  theme_custom() +
  theme(axis.title.x = element_blank()) +
  guides(color_ramp = "none", 
         color = "none")

ggsave("results/mtb-sens.prev.png", width = 12 / cm(1), height = 8 / cm(1))  

# Sensitivity analysis (Outdoor CO2): Comparing 600 ppm vs 400 ppm #

df.sens.co2.sa <- sens.df("South Africa")
df.sens.co2.ch <- sens.df("Switzerland")
df.sens.co2.tz <- sens.df("Tanzania")

df.sens.co2 <- bind_rows(df.sens.co2.sa, df.sens.co2.ch, df.sens.co2.tz) 
  
df.sens.co2 %>%
  mutate(sens = factor(sens, levels = c("South Africa", "Switzerland", "Tanzania"))) %>%
  ggplot(aes(x = sens, y = P, color = country, shape = type)) +
  stat_interval(position = position_dodge(width = 0.7)) +
  stat_summary(geom = "point", fun = "median", 
               position = position_dodge(width = 0.7), size = 2, fill = "white") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), 
                     expand = expansion(add = c(0, 0.0)), 
                     limits = c(0,.5)) +
  scale_color_manual(values = wes_palette("Moonrise2"), guide = "none") +
  scale_shape_manual(values=c(23, 24)) +  # 23 is a square, 24 is a triangle
  labs(y = "Risk of infection (%)", color = "Country", shape = "CO2 Level") +
  theme_custom() +
  theme(axis.title.x = element_blank())

ggsave("results/mtb-sens.co2.png", width = 12 / cm(1), height = 8 / cm(1))  

