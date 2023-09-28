#### Libraries ####

library(tidyverse)
source("utils/settings.r")
source("utils/plotting.r")

#### Data ####

co2 <- readRDS("data-clean/co2.rds")

#### CO2 levels ####

means <- format(round(co2$C.mean), big.mark = ",")

plot.co2 <- co2 %>%
  dplyr::select(country, C) %>%
  unnest() %>%
  mutate(C = ifelse(C > 4000, 4000, C)) %>%
  ggplot(aes(x = C, y = ..density.., fill = country)) +
  facet_wrap(~ country, scales = "free_y") +
  geom_histogram(bins = 30, color = "black") +
  geom_label(data = data.frame(x = 1750, 
                              y = c(7e-04, 6.85e-04, 0.004),
                              country = c("South Africa", "Switzerland", "Tanzania"),
                              label = paste0("Mean=", means, "ppm")),
            mapping = aes(x = x, y = y, label = label), 
            size = 8 / cm(1), hjust = 0, fill = "white", color = "black", label.size = NA) +
  scale_x_continuous(labels = scales::comma, breaks = seq(500, 4000, 1000), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05)), labels = function(x) x * 1000) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(x = expression("CO"[2]*" (ppm)"), y = expression("Density"%*%"10"^-3)) +
  theme_custom() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = 2, size = 10))

plot.co2

ggsave(plot = plot.co2, 
       filename = "results/co2-histogram.png", 
       width = 16,
       height = 8,
       units = "cm")

#### Ventilation ####

# compute
vent <- co2 %>%
  dplyr::select(country, C.daily) %>%
  unnest() %>%
  dplyr::select(country, C.max) %>%
  left_join(data.frame(country = names(vol), vol = vol)) %>%
  left_join(data.frame(country = names(n.class), n = n.class)) %>%
  mutate(Q = ((0.13 * C_a * 1000000) / (C.max - C_o)),
         ACH = ((3600 * Q * n) / vol))

# report
vent %>% 
  group_by(country) %>%
  summarise(across(c(Q, ACH), list(mean = mean, sd = sd))) %>%
  mutate_if(is.numeric, round_k, 2)


