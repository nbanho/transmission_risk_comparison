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
