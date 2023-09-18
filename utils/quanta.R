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
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.ticks =  element_line(),
          legend.text = element_text(size = 8))
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

VER_standing_speaking <- 4.9e-3 
IHR_standing <- c(8.36 + 10.65) / 2 / 1000 * 60 

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
  mutate(p = n / sum(n))

p_activ <- activity$p
names(p_activ) <- activity$activity

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
  if (sum(pa) != 1) {
    stop("Proportions do not add to 1.")
  } 
  
  # if you just sample 1
  if (n == 1) {
    a <- sample(names(pa), pa)
    return(cv(1) * ci[pathogen] * VER[a])
  }
  
  # number of samples per activity
  n_pa <- round(n * pa)
  q <- c(cv(n_pa['breathing']) * ci[pathogen] * VER['sitting, breathing'],
         cv(n_pa['speaking']) * ci[pathogen] * VER['sitting, breathing'],
         cv(n_pa['speaking loudly']) * ci[pathogen] * VER['sitting, speaking loudly'])

  return(q)
  
}

#### Distribution ####

#' for the outbreak scenario we consider the whole distribution
cov_q <- ERq(pathogen = "SARS-CoV-2", n = 1e6, pa = p_activ)
mtb_q <- ERq(pathogen = "Mtb", n = 1e6, pa = p_activ)

data.frame(
  rbind(data.frame(pathogen = "SARS-CoV-2", q = cov_q),
        data.frame(pathogen = "Mtb", q = mtb_q))
) %>%
  ggplot(aes(x = q, y = pathogen, fill = pathogen)) +
  geom_boxplot(alpha = .5) +
  scale_x_log10() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

#' While it is reasonable to consider the whole quanta distribution for an outbreak,
#' it is not plausible to consider the extreme tails when estimating the
#' annual transmission risks. For the annual risk, we are more interested in 
#' the average or median rate. The reasoning is somewhat similar as for the 
#' rebreathed air fraction, which is also computed over the whole time period
#' rather than using a distribution, since it is unlikely that students while 
#' be exposed to CO2 > 2,000ppm for a whole year. Similarly, it is unlikely that 
#' students are exposed to q > 1,000 per hour over a whole year.
#' 
#' Therefore, for the annual scenario we consider the 33%, 50%, and 66% percentile as the
#' low, medium, and high scenario respectively

quantile(cov_q, c(.33, .5, .66))
quantile(mtb_q, c(.33, .5, .66))

#' Note that for Mtb, the low/medium/high scenario are in line with the
#' estimates by Andrews et al considering different infectiousness durations. 


##### Example ####

#' This example shows the modeled annual Mtb transmission risks by country for the
#' scenarios above using, for simplicety, the average CO2 level.

q_ann_mtb <- quantile(mtb_q, c(.33, .5, .66))
scenario <- c("Low", "Medium", "High")
  
df <- tibble(
  country = c("South Africa", "Switzerland", "Tanzania"),
  co2 = c(1626, 1802, 643),
  n = c(30, 20, 50),
  prev_m = c(432, 12, 42),
  prev_l = c(232, 5, 11),
  prev_u = c(632, 20, 73),
  q = c(list(q_ann_mtb), list(q_ann_mtb), list(q_ann_mtb)),
  sc = c(list(scenario), list(scenario), list(scenario)),
  t = 1440
) %>%
  mutate(
    f = (co2 - 400) / 31500,
    prev_s = (prev_u - prev_l) / 2 * qnorm(0.975),
    prev = map2(prev_m, prev_s, function(m, s) rtrunc(1e3, spec = "norm", a = 0, mean = m, sd = s))
  ) %>% 
  unnest(c(q, sc)) %>%
  unnest(c(prev)) %>%
  mutate(I = prev / 100000 * n,
         P = 1 - exp(-f*q*I*t/n))

df %>%
  mutate(country = factor(country, levels = c("South Africa", "Switzerland", "Tanzania")),
         sc = factor(sc, levels = c("Low", "Medium", "High"))) %>%
  ggplot(aes(x = sc, y = P)) +
  stat_interval(aes(color = country, color_ramp = after_stat(rev(.width))), position = position_dodge(width = .5)) +
  stat_summary(aes(x = sc, y = P, color = country), geom = "point", fun = "median", 
               position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
  scale_y_continuous(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0.05)), limits = c(0,1)) +
  scale_x_discrete(labels = c(expression(atop("Low", q*' = 0.8'~'h'^-1)), 
                              expression(atop("Medium", q*' = 3.1'~'h'^-1)),
                              expression(atop("High", q*' = 10.8'~'h'^-1)))) +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  ggdist::scale_color_ramp_continuous() +
  labs(y = "Risk of infection (%)", color_ramp = "Interval", color = "",
       title = "Modeled risk of Mtb transmission",
       subtitle = "Annual risk for t=1,440 school-hours by country",
       caption = "Median as dots, 50%-, 80%-, and 95%-CI as ribbons.") +
  theme_custom() +
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        plot.title.position = "plot",
        legend.box = "vertical") +
  guides(color_ramp = "none", 
         color = guide_legend(order = 1))

ggsave("results/example-mtb-visualization.png", width = 12 / cm(1), height = 8 / cm(1))  
