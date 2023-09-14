#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(truncnorm)
library(LaplacesDemon)
library(tidybayes)

#### Utils ---------------------------------------------------------------------

source("utils/functions-analysis.R")
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

# Main analysis #

f_bar.sa <- f.sa$f_bar
f_bar.ch <- f.ch$f_bar
f_bar.tz <- f.tz$f_bar

# Sensitivity analysis 2 #

f_bar.sens.sa <- f.sa$f_bar.sens
f_bar.sens.ch <- f.ch$f_bar.sens
f_bar.sens.tz <- f.tz$f_bar.sens

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

# Mikzewski standing #

q.logmean_Mikz_sta = log(3.1)
q.logsd_Mikz_sta = (log(430) - log(0.0023)) / (2 * 1.645)
sample.q_Mikz_sta <- rlnorm(n.sample, meanlog = q.logmean_Mikz_sta, sdlog = q.logsd_Mikz_sta)

quant

# Mikzewski standing trucated at 300 # 

sample.q_Mikz_sta.trunc <- rtrunc(n.sample, 
                            spec = "lnorm", 
                            b = 100, 
                            meanlog = q.logmean_Mikz_sta, 
                            sdlog = q.logsd_Mikz_sta)

# Mikzewski resting trucated at 300 # 

q.logmean_Mikz_rest.trunc = log(0.62)
q.logsd_Mikz_rest.trunc = (log(85) - log(0.0045)) / (2 * 1.645)
sample.q_Mikz_rest.trunc <- rtrunc(n.sample, 
                                  spec = "lnorm", 
                                  b = 300, 
                                  meanlog = q.logmean_Mikz_rest.trunc, 
                                  sdlog = q.logsd_Mikz_rest.trunc)

### I (infectious students) ----

# Main analysis #

prev.sa <- rtruncnorm(n.sample, mean = 432, sd = (632 - 232) / (2 * 1.96))
prev.ch <- rtruncnorm(n.sample, a = 0, mean = 12, sd = 3.7)
prev.tz <- rtruncnorm(n.sample, a = 0, mean = 42, sd = (73 - 11) / (2 * 1.96))

I.sa = prev.sa / 100000 * n.sa
I.ch = prev.ch / 100000 * n.ch
I.tz = prev.tz / 100000 * n.tz

# Sensitivity analysis 1 #

prev.sa.sens <- rtruncnorm(n.sample, a = 0, mean = 852, sd = ((1026 - 679) / 2 * 1.96))
prev.ch.sens <- rtruncnorm(n.sample, a = 0, mean = 56.9, sd = 6.4)
prev.tz.sens <- rtruncnorm(n.sample, a = 0, mean = 293, sd = ((358 - 228) / 2 * 1.96))

I.sa.sens = prev.sa.sens / 100000 * n.sa
I.ch.sens = prev.ch.sens / 100000 * n.ch
I.tz.sens = prev.tz.sens / 100000 * n.tz

# Additional analysis #

I.sa.add = 1
I.ch.add = 1 / n.sa * n.ch
I.tz.add = 1 / n.sa * n.tz

#### Results  ------------------------------------------------------------------

data <- tibble(country = c(rep("SA", n.sample), rep("CH", n.sample), rep("TZ", n.sample)),
                I = c(I.sa, I.ch, I.tz),
               I.sens = c(I.sa.sens, I.ch.sens, I.tz.sens),
               I.add = c(rep(I.sa.add, 3000), rep(I.ch.add, 3000), rep(I.tz.add, 3000)),
               n = c(rep(n.sa, n.sample), rep(n.ch, n.sample), rep(n.tz, n.sample)),
               f = c(rep(f_bar.sa, n.sample), rep(f_bar.ch, n.sample), rep(f_bar.tz, n.sample)),
               f.sens = c(rep(f_bar.sens.sa, n.sample), rep(f_bar.sens.ch, n.sample), rep(f_bar.sens.tz, n.sample)),
               q = 0.89,
               t = year,
               t.add = day, 
               P = 1 - exp(-(f * I * q * t) / n),
               P.add = 1 - exp(-(f * I.add * q * t.add) / n),
               P.sens.I = 1 - exp(-(f * I.sens * q * t) / n),
               P.sens.f = 1 - exp(-(f.sens * I * q * t) / n))

# Main analysis #

summary.main <- summarize("P")

ggplot(summary.main, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none") +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 0.3))


ggplot(data, aes(x = country, y = P, fill = country)) +
  geom_violin()

# Additional analysis: fixed I/n #

summary.add <- summarize("P.add")

ggplot(summary.add, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1))

# Sensitivity analysis 1: prevalence in total population #

summary.sens.I <- summarize("P.sens.I")

ggplot(summary.sens.I, aes(x = country, y = middle, fill = country)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "grey") +
  geom_boxplot(aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax), stat = "identity", width = 0.5, alpha = 1, outlier.shape = NA) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = c(0, 0)) +
  scale_fill_manual(values = c("#B0C4DE", "#FFA07A", "#FFDEAD")) +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),  # This ensures the title is centered
        plot.title.position = "panel") +  # This puts the title in the center of the plot area
  coord_cartesian(ylim = c(0, 1))

# Sensitivity analysis 2: Outdoor CO2-level #



