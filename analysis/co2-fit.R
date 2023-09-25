#### Libraries ####
library(tidyverse)
library(fitdistrplus)
library(truncdist)
library(LaplacesDemon)
library(gridExtra)

source(file = "utils/plotting.R")
text_size = 8

#### Data import ####

co2.sa <- readRDS("data-clean/co2-sa.rds") %>% 
  arrange(co2)

co2.ch <- readRDS("data-clean/co2-ch.rds") %>% 
  arrange(co2)

co2.tz <- readRDS("data-clean/co2-tz.rds") %>% 
  arrange(co2)

set.seed(1343)

#### Smoothing the CO2 data ####

# South Africa #

df.sa <- data.frame(
  x = 1:length(co2.sa$co2),
  y = co2.sa$co2)

loess.sa <- data.frame(
  x = 1:length(df.sa$x),
  y = predict(loess(y~x, df.sa, 
                    span = 0.3)),
  method = "loess()"
)

ggplot(loess.sa, aes(x, y)) + 
  geom_point(dat = df.sa, aes(x, y), alpha = 0.2, col = "red") +
  geom_line(col = "blue") +
  facet_wrap(~method) +
  theme_bw(16)

loess.sa %>% #density plot sa smooth
  ggplot(aes(x = y)) +
  geom_density(alpha = .2, kernel = "gaussian", adjust = 4, fill = "darkseagreen3") +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 0.002), expand = c(0,0)) +
  scale_x_continuous(limits = c(400, 3000), expand = c(0,0))  +
  labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
  theme(legend.position = "none") +
  theme_bw()

# Switzerland #

df.ch <- data.frame(
  x = 1:length(co2.ch$co2),
  y = co2.ch$co2)

loess.ch <- data.frame(
  x = 1:length(df.ch$x),
  y = predict(loess(y~x, df.ch, span = 0.1)),
  method = "loess()"
) %>% 
  mutate(school = "Switzerland") #smooting the data using Loess

ggplot(loess.ch, aes(x, y)) + 
  geom_point(dat = df.ch, aes(x, y), alpha = 0.2, col = "red") +
  geom_line(col = "blue") +
  facet_wrap(~method) +
  theme_bw()

# Tanzania #

df.tz <- data.frame(
  x = 1:length(co2.tz$co2),
  y = co2.tz$co2
)

loess.tz <- data.frame(
  x = 1:length(df.tz$x),
  y = predict(loess(y~x, df.tz, span = 0.3)),
  method = "loess()"
)

ggplot(loess.tz, aes(x, y)) + 
  geom_point(dat = df.tz, aes(x, y), alpha = 0.2, col = "red") +
  geom_line(col = "blue") +
  facet_wrap(~method) +
  theme_bw()

loess.tz %>% #density plot tz smooth
  ggplot(aes(x = y)) +
  geom_density(alpha = .2, kernel = "gaussian", adjust = 3.2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 0.0075)) +
  scale_x_continuous(limits = c(400, 3000)) +
  labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
  theme_bw() +
  theme(legend.position = c(0.9,0.9), legend.title = element_blank()) 

#### fitting distribution ####

x <- seq(400, 4000, by = .1) #vector for the plot

# South Africa #

descdist(loess.tz$y, discrete = FALSE)
fit.sa <- fitdistr(loess.tz$y, "lognormal") #get parameters
meanlog.sa <- fit.sa$estimate["meanlog"]
sdlog.sa <- fit.sa$estimate["sdlog"]

distr.sa <- data.frame(co2 = x) %>% 
  mutate(prob = dtrunc(co2, 
                       spec = "lnorm", 
                       a = 400, 
                       b = 4000, 
                       meanlog = meanlog.sa, 
                       sdlog = sdlog.sa),
         country = "South Africa")

# Switzerland #

descdist(loess.ch$y, discrete = FALSE)
fit.ch <- fitdistr(loess.ch$y, "normal")
mean.ch <- fit.ch$estimate["mean"]
sd.ch <- fit.ch$estimate["sd"]

distr.ch <- data.frame(co2 = x) %>% 
  mutate(prob = dtrunc(co2, 
                       spec = "norm", 
                       a = 400, 
                       b = 4000, 
                       mean = mean.ch, 
                       sd = sd.ch),
         country = "Switzerland")

# Tanzania # 

descdist(loess.tz$y, discrete = FALSE, boot = 100) 
fit.tz <- fitdistr(loess.tz$y, "t") #get parameters
location.tz <- fit.tz$estimate["m"]
scale.tz <- fit.tz$estimate["s"]
df_par.tz <- fit.tz$estimate["df"]

distr.tz <- data.frame(co2 = x) %>% 
  mutate(prob = dtrunc(co2, 
                       spec = "st", 
                       a = 400, 
                       b = 4000, 
                       mu = location.tz, 
                       sigma = scale.tz,
                       nu = df_par.tz),
         country = "Tanzania")

#### Plotting ####

## i adjusted the parameters for South Africa and Tanzania slightly ##

# South Africa # 

plot_co2.sa <-  ggplot(data = loess.sa, aes(x = y)) +
  geom_histogram(binwidth = 100, color = "black", fill = "#B0C4DE", aes(y = ..density..)) +
  scale_y_continuous(limits = c(0, 0.001), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
  coord_cartesian(xlim = c(400,3500))+
  theme_bw2()+
  stat_function(fun = function(x) dtrunc(x, spec = "lnorm", a = 400, b = 4000, meanlog = 7.25, sdlog = 0.5), linewidth = 1) +
  ggtitle("South Africa") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2840, y = 0.000965, label = expression(paste("Log-normal (", mu, " = 7.3, ", sigma, " = 0.5)")), size = 3, hjust = 0.7, vjust = 1)

print(plot_co2.sa)

# Switzerland #

plot_co2.ch <- loess.ch %>% 
  ggplot(aes(x = y)) +
  geom_histogram(binwidth = 100, color = "black", fill = "#FFA07A",  aes(y = ..density..)) +
  scale_y_continuous(limits = c(0, 0.001), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  labs(x = expression(CO[2]*" (ppm)"), y = "") +
  coord_cartesian(xlim = c(400,3500))+
  theme_bw2()+ 
  stat_function(fun = function(x) dtrunc(x, spec = "norm", a = 400, b = 4000, mean = 1802, sd = 546), linewidth = 1) +
  ggtitle("Switzerland") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2840, y = 0.00095, label = expression(paste("Normal (", mu, " = 1802, ", sigma, " = 546)")), size = 3, hjust = 0.7, vjust = 1)

plot_co2.ch

# Tanzania #

plot_co2.tz <- loess.tz %>% 
  ggplot(aes(x = y)) +
  geom_histogram(binwidth = 100, color = "black", fill = "#FFDEAD", aes(y = after_stat(density))) +
  scale_y_continuous(limits = c(0, 0.006), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
  coord_cartesian(xlim = c(400,3500))+
  stat_function(fun = function(x) dtrunc(x, spec = "st", a = 400, b = 4000, mu = 594, sigma = 80, nu = 1.55), size = 1) +
  ggtitle("Tanzania") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2850, y = 0.0058, label = expression(paste("Student-t (", mu, " = 594, ", sigma, " = 80, ", nu, " = 1.55)")), size = 3, hjust = 0.75, vjust = 1)+
  theme_bw2() 

print(plot_co2.tz)

# together #

plot.co2 <- grid.arrange(plot_co2.sa,plot_co2.ch,plot_co2.tz, ncol = 2)

# save #

ggsave(plot = plot.co2, 
       filename = "results/co2/co2-distr.png", 
       width = 8 ,
       height = 5 
       )

library(wesanderson)

## Plotting raw data only ##
mean(co2.sa$co2)

plot.sa <- co2.sa %>% 
  ggplot(aes(x = co2)) +
  geom_histogram(aes(y =..density..), color = "black", fill = wes_palette("Moonrise2")[1], bins = 30) +
  scale_y_continuous(limits = c(0, 0.001), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  coord_cartesian(xlim = c(400,3500))+
  labs(x = "", y = "Density") +
  theme_classic() +
  ggtitle("South Africa") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2840, y = 0.000965, label = "Mean = 1,626 ppm", size = 3, hjust = 0.7, vjust = 1)

mean(co2.ch$co2)

plot.ch <- co2.ch %>% 
  ggplot(aes(x = co2)) +
  geom_histogram(aes(y =..density..), color = "black", fill = wes_palette("Moonrise2")[2], bins = 30) +
  scale_y_continuous(limits = c(0, 0.001), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  labs(x = expression(CO[2]*" (ppm)"), y = "") +
  coord_cartesian(xlim = c(400,3500))+
  theme_classic() +
  ggtitle("Switzerland") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2840, y = 0.000965, label = "Mean = 1,802 ppm", size = 3, hjust = 0.7, vjust = 1)

mean(co2.tz$co2)

plot.tz <- co2.tz %>% 
  ggplot(aes(x = co2)) +
  geom_histogram(aes(y =..density..), color = "black", fill = wes_palette("Moonrise2")[3], bins = 30) +
  scale_y_continuous(limits = c(0, 0.006), expand = c(0,0), labels = function(x) sprintf("%.4f", x)) +
  scale_x_continuous(limits = c(400, 4000), expand = c(0,0), labels = scales::comma) +
  labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
  coord_cartesian(xlim = c(400,3500))+
  theme_classic() +
  ggtitle("Tanzania") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2850, y = 0.0058, label = "Mean = 648 ppm", size = 3, hjust = 0.75, vjust = 1)

# together #

plot.co2.raw <- grid.arrange(plot.sa,plot.ch,plot.tz, ncol = 2)

# save #

ggsave(plot = plot.co2.raw, 
       filename = "results/co2/co2-distr_raw.png", 
       width = 8 ,
       height = 5 
)

