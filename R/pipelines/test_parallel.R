options(scipen = 9999999)
options(digits = 20)
library("tidyverse")
library("stringr")
library("foreach")
library("doParallel")
source("R/data_processing/calculate_cosine_similarity.R")

# Set the number of cores to use for parallel processing
num_cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
registerDoParallel(cores = num_cores / 2)

habitats <- dir("data/performance_fastspar_iterations/", full.names = TRUE)
data_frames <- list()

# function to convert the data frame to numeric avoiding errors
convert_to_numeric <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.numeric)
  df
}

foreach(habitat_path = habitats, .packages = c("tidyverse", "stringr")) %dopar% {


  habitat_name <- basename(habitat_path)

  # intialized the data frame of the habitat results
  habitat_data_frame <- data.frame()

  # order the iterations ascendingly
  iterations <- dir(habitat_path, full.names = TRUE)

  orderer <- as.numeric(str_extract(iterations, "\\d+"))
  iterations <- iterations[order(orderer)]

  # get iteration folders
  foreach(iteration_path = iterations, .packages = c("tidyverse", "stringr")) %dopar% {

    iteration_name <- basename(iteration_path)
    print(paste0("checking results for:", habitat_name, "iteration:", iteration_name))

    # get the tables
    tables <- dir(iteration_path, pattern = "cor_", full.names = TRUE)

    orderer <- as.numeric(str_extract(tables, "\\d+"))
    tables <- tables[order(orderer)]

    matrices <- list()

    foreach(table_path = tables, .packages = c("tidyverse")) %dopar% {

      table_name <- basename(table_path)
      print(paste0("table:", table_name))
      table_name <- basename(table_path)

      # treat from df to matrix
      table <- as.matrix(read_tsv(table_path, show_col_types = FALSE))
      table <- as.data.frame(table)
      table[1] <- NULL
      colnames(table) <- NULL
      table <- convert_to_numeric(table)
      table <- as.matrix(table)

      # store the matrix
      matrices[[table_name]] <- table

    }

    results_similarities <- compute_matrices_similarity(matrices)
    results_similarities <- as.vector(results_similarities[upper.tri(results_similarities, diag = FALSE)])
    results_df_parcial <- data.frame(iteration_name = results_similarities)
    colnames(results_df_parcial) <- c(iteration_name)

    habitat_data_frame <- bind_cols(habitat_data_frame, results_df_parcial)

  }

  result_path <- "data/summaries/performance_fastspar_iterations/"
  if (!file.exists(result_path)) {
    dir.create(result_path)
  }
  write.csv(habitat_data_frame, paste0(result_path, habitat_name, ".csv"))

  data_frames[[habitat_name]] <- habitat_data_frame
}
