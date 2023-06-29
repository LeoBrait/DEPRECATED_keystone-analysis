options(scipen = 9999999)
options(digits = 20)
library("tidyverse")
library("stringr")
source("R/data_processing/calculate_cosine_similarity.R")

args <- commandArgs(trailingOnly = TRUE)
frame_analysis <- as.character(args[1])

habitats <- dir(paste0(
  "data/", frame_analysis, "performance_fastspar_iterations/",
  full.names = TRUE))

data_frames <- list()
results_df <- data.frame()

#function to convert the data frame to numeric avoiding errors
convert_to_numeric <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.numeric)
  df
}


for (habitat_path in habitats){


  habitat_name <- basename(habitat_path)

  if(file.exists(
    paste0(
      "data/",
         frame_analysis,
            "summaries/performance_fastspar_iterations/",
               habitat_name,
                 ".csv"
                 ))){
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
  result_path <- paste0(
    "data/",
      frame_analysis,
        "summaries/performance_fastspar_iterations/")
  if (!file.exists(result_path)) {
    dir.create(result_path)}
  
  write.csv(
    results_df,
    paste0(result_path, habitat_name, ".csv"),
    row.names = FALSE)
  
  data_frames[[habitat_name]] <- results_df
}}




