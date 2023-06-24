options(scipen = 9999999)
library("doParallel")
library("tidyverse")
library("stringr")
library("foreach")
library("readr")

source("R/data_processing/calculate_cosine_similarity.R")

# Manual test
habitat <- dir("data/performance_fastspar_iterations/", full.names = TRUE)


iter_names <- dir("data/performance_fastspar_iterations/animal_host-associated.animal_feces/")
numeric <- as.numeric(iter_names) %>% sort()
iter_names <- as.character(numeric)


# Set up parallel processing
cores <- detectCores() - 10
cl <- makeCluster(cores)
registerDoParallel(cl)

# Function to compute similarity matrix for a given habitat and iteration
compute_similarity <- function(habitat_path, iteration_name) {

  full_path <- paste0(habitat_path, "/", iteration_name)
  iter_sets_paths <- dir(full_path, pattern = ".cor", full.names = TRUE)
  
  # Load the tables
  iter_sets <- list()
  tables <- list()
  
  for (i in 1:length(iter_sets_paths)) {
    tables[[i]] <- as.matrix(readr::read_tsv(iter_sets_paths[i]))
    iter_sets[[i]] <- as.data.frame(tables[[i]])
    colnames(iter_sets[[i]]) <- NULL
    iter_sets[[i]][1] <- NULL
  }
  
  convert_to_numeric <- function(df) {
    char_cols <- sapply(df, is.character)
    df[char_cols] <- lapply(df[char_cols], as.numeric)
    df
  }
  
  # Convert data frames to numerical matrices
  matrix_list <- lapply(iter_sets, function(df) as.matrix(convert_to_numeric(df)))
  
  # Compute similarity matrix
  results_similarities <- compute_matrices_similarity(matrix_list)
  
  return(as.vector(results_similarities[upper.tri(results_similarities, diag = FALSE)]))
}

# Compute similarities for each habitat and iteration using parallel processing
results <- foreach(path = habitat, .combine = rbind) %:%
           foreach(iteration = iter_names, .combine = rbind) %dopar% {
             compute_similarity(path, iteration)
           }

# Stop parallel processing and clean up
stopCluster(cl)
registerDoSEQ()

# Create a dataframe to store the results
results_df <- as.data.frame(results)
colnames(results_df) <- iter_names

# Print the dataframe
print(results_df)


###foreach
