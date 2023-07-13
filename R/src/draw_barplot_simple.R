#' @title draw_barplot_simple
#' @description Draw a barplot with error bars.
#' @param data the data frame with the data to plot.
#' @param title the title of the plot.
#' @param x_var the variable to plot in the x axis.
#' @param y_var the variable to plot in the y axis.
#' @param error_var the variable to plot as error bars.
#' @param title_y the title of the y axis.
#' @param title_x the title of the x axis.
#' @param legend_title the title of the legend.
#' @param legend_position the position of the legend.
#' @param breaks the breaks of the x axis.
#' @param break_labels the labels of the breaks in the y axis.
#' @param colors the colors of the bars.
draw_barplot_simple <- function(
  data,
  title = "",
  x_var,
  y_var,
  error_var = NULL,
  title_y = "",
  title_x = "",
  legend_title = "",
  legend_position = "none",
  breaks = NULL,
  break_labels = NULL,
  colors) {


  barplot <- ggplot(data) +

      ## Theme and legend
      theme_pubr() +
      theme(text = element_text(size = unit(9, "cm")),
      strip.background = element_blank(),
      legend.position = legend_position,
      legend.text = element_text(size = unit(9, "cm")),
      legend.title = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(face = "bold", size = unit(9, "cm")),
      axis.title.y = element_text(face = "bold", , size = unit(9, "cm"))) +
      guides(fill = guide_legend(
          title = legend_title, size = unit(9, "cm"), shape = 16, nrow = 2)) +

      ## Error bars
      geom_errorbar(
        aes(
          x = .data[[x_var]],
          ymin = .data[[y_var]] - .data[[error_var]],
          ymax = .data[[y_var]] + .data[[error_var]]
          ),
        width = 0.3,
        colour = "black",
        alpha = 0.9,
        linewidth = 0.5) +

      ## Bars
      geom_bar(
        aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]]),
        stat = "identity",
        width = .5) +

      ## legend
      scale_fill_manual(values = colors) +
      theme(
        legend.title = element_text(
          angle = 0,
          size = unit(9, "cm"),
          face = "bold"),
      legend.spacing.x = unit(0.3, "cm"),
      legend.spacing.y = unit(0.3, "cm")) +
      
      ## Titles
      ggtitle(title) +
      theme(
        plot.title = element_text(
          size = unit(9, "cm"),
          face = "bold",
          hjust = 0.5)) +
      labs(x = title_x, y = title_y) +

      ##Breaks
      scale_y_continuous(breaks = breaks, labels = break_labels)

return(barplot)
}