#### Plotting ####

library(wesanderson)

text_size = 8
update_geom_defaults("text", list(size = text_size))

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

#### Risk plots ####

plot_risk <- function(df) {
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
    stat_summary(data = df,
                 mapping = aes(x = scenario, y = P, color = country), geom = "point", fun = "median", 
                 position = position_dodge2(width = .5), size = 2, shape = 23, fill = "white") +
    scale_y_sqrt(labels = scales::percent_format(suffix = ""), expand = expansion(add = c(0, 0)), 
                       limits = c(0,1), breaks = seq(0, 1, .2)^2) +
    scale_x_discrete(labels = c("Low", "Medium", "High"), expand = expansion(add = c(.33, .33))) +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2")) +
    labs(y = "Monthly risk of infection (%, square-root scale)", 
         x = "Activity level in the classroom",
         color = '',
         fill = '') +
    theme_custom() +
    theme(legend.position = "top",
          plot.title.position = "plot",
          legend.box = "vertical") 
}




plot_risk.sens <- function(df) {
  df %>% 
    ggplot(mapping = aes(x = scenario, y = P, fill = country, shape = type)) +
    geom_errorbar(data = df %>% 
                    dplyr::select(country, scenario, P, type) %>% 
                    group_by(country, scenario, type) %>% 
                    median_qi(),
                  mapping = aes(ymin = .lower, ymax = .upper, color = country),
                  position = position_dodge(width = .5),
                  width = .3) +
    geom_boxplot(data = df, 
                 position = position_dodge(width = .5),
                 outlier.shape = NA, coef = 0, width = 0.3) +
    stat_summary(data = df,
                 mapping = aes(x = scenario, y = P, group = country), geom = "point", fun = "median", 
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
    labs(y = "Monthly risk of infection (%, square-root scale)", 
         x = "",
         color = '',
         fill = '') +
    theme_custom() +
    theme(legend.position = "top",
          plot.title.position = "plot",
          legend.box = "vertical") 
}


