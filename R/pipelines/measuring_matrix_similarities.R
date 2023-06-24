options(scipen = 9999999)
library("tidyverse")
library("stringr")
source("R/data_processing/calculate_cosine_similarity.R")

#manual test
habitat <- dir("data/performance_fastspar_iterations/animal_host-associated.habitat", full.names = TRUE)
orderer <- as.numeric(str_extract(habitat, "\\d+"))
habitat <- habitat[order(orderer)]


iter_names <- dir("data/performance_fastspar_iterations/animal_host-associated.habitat")
iter_names <- iter_names[order(orderer)]


for (i in 1:length(habitat)) {

  iter_sets_paths <- list()
  
  for(j in 1:length(iter_names)) {
    iter_sets_paths[[iter_names[j]]] <- dir(habitat[j], pattern = ".cor", full.names = TRUE)
  }

}

#load the tables
iter_sets <- list()
tables <-  list()
for (i in 1:length(iter_sets_paths)) {
  
  iter_sets[[ iter_names[i] ]] <- list()
  
  for (j in 1:length(iter_sets_paths[[i]])) {
    
    tables[[j]] <- as.matrix(read_tsv(iter_sets_paths[[i]][j]))
    iter_sets[[ iter_names[i] ]] [[j]] <- as.data.frame(tables[[j]])
    colnames(iter_sets[[ iter_names[i] ]] [[j]]) <- NULL
    iter_sets[[ iter_names[i] ]] [[j]][1] <- NULL

  }}

convert_to_numeric <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.numeric)
  df
}


# Convert data frames to numerical matrices
matrix_list <- lapply(iter_sets$"16000", function(df) as.matrix(convert_to_numeric(df)))

results_simiralities <- compute_matrices_similarity(matrix_list)

results_simiralities <- as.vector(results_simiralities[upper.tri(results_simiralities, diag = FALSE)])
