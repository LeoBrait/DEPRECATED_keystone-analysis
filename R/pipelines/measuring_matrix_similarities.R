options(scipen = 9999999)
options(digits = 20)
library("tidyverse")
library("stringr")
library("ggpubr")
library("cowplot")
source("R/data_processing/calculate_cosine_similarity.R")

#function to convert the data frame to numeric avoiding errors
convert_to_numeric <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.numeric)
  df
}


args <- commandArgs(trailingOnly = TRUE)
frame_analysis <- as.character(args[1])
summary_path <- paste0(
    "data/",
      frame_analysis,
        "/summaries/performance_fastspar_iterations/")

plot_path <- paste0(
    "results/",
      frame_analysis,
        "/performance_fastspar_iterations/")

if(!file.exists("results/")){
  dir.create("results/")}

if (!file.exists(summary_path)) {
    dir.create(summary_path)}

if (!file.exists(plot_path)) {
    dir.create(plot_path)}


habitats <- dir(paste0(
  "data/", frame_analysis, "/performance_fastspar_iterations/"),
  full.names = TRUE)


data_frames <- list()
results_df <- data.frame()





for (habitat_path in habitats){
    
  habitat_name <- basename(habitat_path)
  print(paste0("processing: ", habitat_name))

  if(file.exists(
    paste0(summary_path, habitat_name,".csv"))){
    print(paste0("skipping: ", habitat_name))

  }else{

  #intialized the data frame of the habitat results
  habitat_data_frame <- data.frame()

  #order the iterations ascendingly
  iterations <- dir(habitat_path, full.names = TRUE)

  orderer <- as.numeric(str_extract(iterations, "\\d+"))
  iterations <- iterations[order(orderer)]

  #get iteration folders
  for(iteration_path in iterations){

  iteration_name <- basename(iteration_path)
  print(paste0(
    "checking results for: ", habitat_name, "iteration: ", iteration_name))
    
  #get the tables
  tables <- dir(iteration_path, pattern = "cor_", full.names = TRUE)

  orderer <- as.numeric(str_extract(tables, "\\d+"))
  tables <- tables[order(orderer)]

  matrices <- list()


    for(table_path in tables){

      table_name <- basename(table_path)
      print(paste0("table: ", table_name))
      table_name <- basename(table_path)

      #treat from df to matrix
      table <- as.matrix(read_tsv(table_path, show_col_types = FALSE))
      table <- as.data.frame(table)
      table[1] <- NULL
      colnames(table) <- NULL
      table <- convert_to_numeric(table)
      table <- as.matrix(table)

      #store the matrix
      matrices[[table_name]] <- table

    }
  
  results_similarities <- compute_matrices_similarity(matrices)
  results_similarities <- as.vector(results_similarities[upper.tri(results_similarities, diag = FALSE)])
  results_df_parcial <- data.frame(iteration_name = results_similarities)
  colnames(results_df_parcial) <- c(iteration_name)

  #This prevents an error while binding the data frames
  if (nrow(results_df) == 0) {
       results_df <- results_df_parcial
   } else {
       results_df <- cbind(results_df, results_df_parcial)
   }

  }

  
  write.csv(
    results_df,
    paste0(summary_path, habitat_name, ".csv"),
    row.names = FALSE)
  

  results_df <- data.frame()
  
   }
}

plot_list <- list()
#habitats <- habitats[1]
for (habitat_path in habitats){

  habitat_name <- basename(habitat_path)

  data <- read_csv(paste0(summary_path, habitat_name, ".csv"))

  data_long <- tidyr::gather(data, Iterations, similarity)

  habitat_name <- str_replace_all(habitat_name, "_", " ")
  habitat_name <- str_replace_all(habitat_name, "\\.", " ")


  plot_list[[habitat_name]] <- ggplot(
    data_long, aes(x = Iterations, y = similarity, group = 1)) +
    theme_pubr() +
      theme(text = element_text(size = unit(9, "cm")),
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(face = "bold", size = unit(9, "cm")),
      axis.title.y = element_text(face = "bold", , size = unit(9, "cm"))) +
      theme(text = element_text(family = "arial")) +
      geom_line() +
      labs(x = "Iterations", y = "Simmilarity (%)") +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 0.95)) +
      
      ## Titles
      ggtitle(habitat_name) +
      theme(
        plot.title = element_text(
          size = unit(9, "cm"),
          face = "bold",
          hjust = 0.5)) +

      scale_x_discrete(limits = colnames(data)) +
      scale_y_continuous(n.breaks = 2, labels =c("99", "100"))
  

}

panel <- plot_grid(plotlist = plot_list, ncol = 3)
plot(panel)

ggsave(
  paste0(plot_path, "panel.png"),
  plot = panel,
  width = 19,
  height = 10,
  units = "cm",
  dpi = 300)

