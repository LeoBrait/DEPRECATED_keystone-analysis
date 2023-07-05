
library("ggplot2")
library("ggpubr")

  draw_lineplot <- function(
    data,
    x_var,
    y_var,
    habitat_name = NULL,
    x_title = NULL,
    y_title = NULL){

    plot <- ggplot(
        data = data,
        aes(x = data[[x_var]], y = data[[y_var]], group = 1)) +
        theme_pubr() +
        theme(
            text = element_text(size = unit(9, "cm")),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(face = "bold", size = unit(9, "cm")),
            axis.title.y = element_text(
                face = "bold", , size = unit(9, "cm"))) +
        theme(text = element_text(family = "Arial")) +
        geom_line() +
        labs(x = x_title, y = y_title) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 0.95)) +

        ## Titles
        ggtitle(habitat_name) +
        theme(plot.title = element_text(size = unit(9, "cm"), face = "bold")) +
        scale_y_continuous(
        breaks = c(min(data[[y_var]]), max(data[[y_var]])),
        labels = c(
          floor(min(data[[y_var]]) * 100), "100")
        )
    return(plot)
  }

