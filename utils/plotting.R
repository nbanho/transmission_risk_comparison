#### Plotting ####

plot_co2 <- function(data_frame) {
  ggplot(data_frame, aes(x = y)) +
    geom_histogram(binwidth = 100, color = "black", fill = "#FFDEAD", aes(y = after_stat(density))) +
    scale_x_continuous(limits = c(400, 4000), expand = c(0, 0), labels = scales::comma) +
    labs(x = expression(CO[2]*" (ppm)"), y = "Density") +
    coord_cartesian(xlim = c(400, 3500)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
}

theme_bw2 <- function () { 
  theme_bw(base_size = text_size, base_family = "sans") %+replace% 
    theme(
      axis.text = element_text(size = text_size),
      axis.title = element_text(size = text_size),
      plot.title = element_text(size = text_size + 2, face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 5, 0))
    )
}

cm <- function(x) {
  x / 2.54
}

#### Results -------------------------------------------------------------------

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
               q = q.med,
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


