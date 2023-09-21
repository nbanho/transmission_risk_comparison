#### Libraries ####

library(tidyverse)
library(reshape2)
library(LaplacesDemon)
library(tidybayes)

# for plotting
library(wesanderson)
text_size = 8
update_geom_defaults("text", list(size = text_size))
theme_custom <- function() {
  theme_minimal() %+replace% 
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size),
          axis.title = element_text(size = text_size),
          plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 5, 0)),
          strip.text = element_text(size = text_size),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.ticks =  element_line(),
          legend.text = element_text(size = 8))
}


# for rounding so that the vector sums to 1
smart.round <- function(x) {
  if (isTRUE(all.equal(sum(x), 1))) {
    x <- x * 100
  }
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}


#### Derivation ####

#' References:
#' - Mikszewski: Mikszewski et al. (2022) in Geoscience frontiers (doi.org/10.1016/j.gsf.2021.101285)
#' - Adams: Adams (1993) [Report] Measurement of Breathing Rate and Volume in Routinely Performed Daily Activities.
#' - Stadnytskyi: Stadnytskyi et al. (2020) in PNAS (doi.org/10.1073/pnas.2006874117)
#' - Morawska: Morawska et al. (2009) in Journal of Aerosol Science (doi.org/10.1016/j.jaerosci.2008.11.002)
#' - Buonanno: Buonanno et al. (2020) in Environment International (doi.org/10.1016/j.envint.2020.106112)
#' - Andrews: Andrews et al. (2014) in JID (doi.org/10.1093/infdis/jiu138)

#' Mikszewski estimates the quanta generation rate distribution using 
#' the inhalation rates (IR) for resting, standing, and light activity (slow walking) 
#' from Adams (1993). 
#' However, students are mostly sitting, which has a lower IR. Therefore, we want to
#' re-estiamte the quanta distribution using the IR for sitting. Note that we still assume,
#' like Adams, that adolescents having similar IRs as adults, i.e. we use the 
#' average IR of adults female/male. 

IR <- c(c(7.12 + 8.98) / 2 / 1000 * 60,
        c(7.72 + 9.30) / 2 / 1000 * 60,
        c(8.36 + 10.65) / 2 / 1000 * 60,
        c(20.32 + 24.13) / 2 / 1000 * 60)
names(IR) <- c("lying (resting)", "sitting", "standing", "walking slowly (light activity)")

#' Mikszewski computes the droplet volume emission rate VER as the product
#' of the IR and the droplet volume concentration Vd, using the droplet 
#' emission rates from Stadnytskyi. Stadnytskyi estimates the transmission rates
#' only for speaking, and (I assume) Mikszewski scales them to other respiratory 
#' activities (breathing, loud speaking), using the multipliers derived 
#' from the estimates for Vd from Morawska as done in Buoanno. 
#' 
#' Therefore, we re-scale the VER_M reported by Mikszewski by multiplying them with 
#' IR_sitting / IR_ALM, where IR_ALM is the IR of the activity level that was 
#' actually used by Mikszewski.  

VER_M <- c(9.8e-4, 4.9e-3, 8.3e-2)
names(VER_M) <- c("resting, oral breathing", "standing, speaking", "light activity, speaking loudly")

VER <- c(IR["sitting"] / IR["lying (resting)"] * VER_M["resting, oral breathing"],
         IR["sitting"] / IR["standing"] * VER_M["standing, speaking"],
         IR["sitting"] / IR["walking slowly (light activity)"] * VER_M["light activity, speaking loudly"])
names(VER) <- c("sitting, breathing", "sitting, speaking", "sitting, speaking loudly")

#' To compute the quanta rate, we use the vira the viral load input data 
#' and conversion factors reported by Mikszewski in Table 1

cv_cov <- function(n) rlnorm(n, log(10^5.6), log(10^1.2)) # RNA copies mL-1
cv_mtb <- function(n) rlnorm(n, log(10^5.5), log(10^1.3)) # CFU mL-1
ci <- c(1.4e-3, 2.0e-3)
names(ci) <- c("SARS-CoV-2", "Mtb")

#' Finally, we use the data from our Swiss study to inform the average prop. of activity levels
#' in the classroom. In the study, we tracked activities at 10min intervals and
#' categorized them into "silent working" (breathing), "quiet working, speaking"
#' (speaking), "loud speaking" (loud speaking).   

activity <- read.csv("data-raw/activity.csv") %>%
  mutate(class = ifelse(grepl("3a", record_id), "3a", "3b")) %>%
  dplyr::select(class, date, matches("morning|midday|afternoon")) %>%
  dplyr::select(class, date, matches("___1"), matches("___2"), matches("___3")) %>%
  melt(c("class", "date")) %>%
  mutate(activity = ifelse(grepl("___1", variable), "breathing", 
                           ifelse(grepl("___2", variable), "speaking", 
                                  "speaking loudly"))) %>% 
  mutate(time = stringi::stri_extract(variable, regex = "\\d{3,4}")) %>%
  mutate(value = ifelse(value == 1, T, F)) %>%
  dplyr::select(class, date, time, activity, value) %>%
  group_by(class, date, time) %>%
  filter(any(value)) %>% # exclude times without lesson in the classroom
  ungroup() %>%
  group_by(activity) %>%
  dplyr::summarize(n = sum(value)) %>%
  ungroup() %>%
  mutate(p = n / sum(n)) %>%
  mutate(activity = factor(activity, levels = c("breathing", "speaking", "speaking loudly"))) %>%
  arrange(activity)

p_activ <- activity$p
names(p_activ) <- activity$activity

round(p_activ * 100)

#' No we can re-estimate the quanta emission rates depending on activity levels
#' using the predictive approach by Buonanno (see Eq. (1) in Mikszewski):

#' Quanta emission rate sampling distributions ERq for students sitting in the classroom 
#' depending on activity level
#' 
#' @param pathogen character for the pathogen name
#' @param n scalar for the total number of samples
#' @param pa named numeric vector ('breathing', 'speaking', 'speaking loudly') 
#' 
#' @details If n==1, then the activity level will be sampled according to pa. 

ERq <- function(pathogen = "SARS-CoV-2", n, pa = p_activ) {
  
  # use pathogen-specific viral load distribution
  if (pathogen == "SARS-CoV-2") {
    cv <- cv_cov
  } else if (pathogen == "Mtb") {
    cv <- cv_mtb
  } else {
    stop("Pathogen undefined.")
  }
  
  # check activity level names
  sit_a <- paste0("sitting, ", c('breathing', 'speaking', 'speaking loudly'))
  sit_names_pa <- paste0("sitting, ", names(pa))
  if (!all(sit_names_pa %in% sit_a)) {
    stop("Activity levels undefined")
  }
  
  # check that proportions sum to 1
  if (!isTRUE(all.equal(sum(pa), 1))) {
    stop("Proportions do not add to 1.")
  } 
  
  # if you just sample 1
  if (n == 1) {
    a <- sample(names(pa), pa)
    return(cv(1) * ci[pathogen] * VER[a])
  }
  
  # number of samples per activity
  n_pa <- smart.round(n * pa)
  q <- c(cv(n_pa['breathing']) * ci[pathogen] * VER['sitting, breathing'],
         cv(n_pa['speaking']) * ci[pathogen] * VER['sitting, speaking'],
         cv(n_pa['speaking loudly']) * ci[pathogen] * VER['sitting, speaking loudly'])

  return(q)
  
}

#### Distribution ####

n.sample <- 1e6

#' Summaries by activity
set.seed(1)

p_breath <- c(1, 0, 0)
p_speak <- c(0, 1, 0)
p_loud_speak <- c(0, 0, 1)
names(p_breath) <- names(p_speak) <- names(p_loud_speak) <- names(p_activ)

round(quantile(ERq(n = n.sample, pa = p_breath), c(0.5, 0.25, 0.75)), 1)
round(quantile(ERq(n = n.sample, pa = p_speak), c(0.5, 0.25, 0.75)), 1)
round(quantile(ERq(n = n.sample, pa = p_loud_speak), c(0.5, 0.25, 0.75)), 1)

round(quantile(ERq("Mtb", n = n.sample, pa = p_breath), c(0.5, 0.25, 0.75)), 1)
round(quantile(ERq("Mtb", n = n.sample, pa = p_speak), c(0.5, 0.25, 0.75)), 1)
round(quantile(ERq("Mtb", n = n.sample, pa = p_loud_speak), c(0.5, 0.25, 0.75)), 1)

#' Weighted activity summary

p_low <- c("breathing" = 0.7, "speaking" = 0.25, "speaking loudly" = 0.05)
p_med <- c("breathing" = 0.5, "speaking" = 0.4, "speaking loudly" = 0.1)
p_high <- c("breathing" = 0.3, "speaking" = 0.5, "speaking loudly" = 0.2)

q_cov_low <- ERq(n = n.sample, pa = p_low)
q_cov_med <- ERq(n = n.sample, pa = p_med)
q_cov_high <- ERq(n = n.sample, pa = p_high)

round(quantile(q_cov_low, c(.5, .25, .75)), 1)
round(quantile(q_cov_med, c(.5, .25, .75)), 1)
round(quantile(q_cov_high, c(.5, .25, .75)), 1)

q_mtb_low <- ERq("Mtb", n.sample, pa = p_low)
q_mtb_med <- ERq("Mtb", n.sample, pa = p_med)
q_mtb_high <- ERq("Mtb", n.sample, pa = p_high)

round(quantile(q_mtb_low, c(.5, .25, .75)), 1)
round(quantile(q_mtb_med, c(.5, .25, .75)), 1)
round(quantile(q_mtb_high, c(.5, .25, .75)), 1)


##### Mtb Example ####

#' This example shows the modeled annual Mtb transmission risks by country.
#' 
#' Further assumptions:
Ca <- 31500
Co <- 400
t_ann <- 919

df <- tibble(
  country = c("South Africa", "Switzerland", "Tanzania"),
  n = c(30, 20, 50),
  C = c(1626, 1802, 643),
  prev_m = c(432, 12, 42),
  prev_l = c(232, 5, 11),
  prev_u = c(632, 20, 73),
  f = (C - Co) / Ca,
  prev_s = (prev_u - prev_l) / (2 * qnorm(0.975)),
  q =  rep(list(c(q_mtb_low, q_mtb_med, q_mtb_high)), each = 3),
  scenario = rep(list(c(rep("Low", n.sample), rep("Medium", n.sample), rep("High", n.sample)))),
  t = t_ann
) 

rI <- function(n, pm, ps, nstud) {
  rtrunc(n, spec = "norm", a = 0, mean = pm, sd = ps) / 1e5 * nstud
}

df <- df %>%
  mutate(I = pmap(list(prev_m, prev_s, n), function(m, s, n) rep(rI(n.sample, m, s, n), 3))) %>%
  unnest(c(I, q, scenario)) %>%
  mutate(P = 1 - exp(-f*q*I*t/n)) %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania")),
         scenario = factor(scenario, levels = c("Low", "Medium", "High")))

ggplot(mapping = aes(x = scenario, y = P, fill = country)) +
  geom_errorbar(data = df %>%
                  dplyr::select(country, scenario, P) %>%
                  group_by(country, scenario) %>%
                  median_qi(),
                mapping = aes(ymin = .lower, ymax = .upper, color = country),
                position = position_dodge(width = .5),
                width = .3) +
  geom_boxplot(data = df, 
               position = position_dodge(width = .5),
               outlier.shape = NA, coef = 0, width = 0.3) +
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

ggsave("results/example-mtb-visualization.png", width = 12 / cm(1), height = 8 / cm(1))  
