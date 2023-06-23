sci.pe
library("tidyverse")
source("R/data_processing/calculate_cosine_similarity.R")

#manual test
aqueous_humour <- dir("data/performance_fastspar_iterations/animal_host-associated.aqueous_humour", full.names = TRUE)
iter_names <- dir("data/performance_fastspar_iterations/animal_host-associated.aqueous_humour")

for (i in 1:length(aqueous_humour)) {

  iter_sets_paths <- list()
  
  for(j in 1:length(iter_names)) {
    iter_sets_paths[[iter_names[j]]] <- dir(aqueous_humour[j], pattern = ".cor", full.names = TRUE)
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
matrix_list <- lapply(iter_sets$"300", function(df) as.matrix(convert_to_numeric(df)))

x <- as.matrix(matrix_list[[1]])


x <- cosine_similarity(matrix_list[[1]], matrix_list[[2]])


flattened1 <- c(matrix_list[[1]])
flattened2 <- c(matrix_list[[2]])
dot_product <- sum(flattened1 * flattened2)
norm_product <- sqrt(sum(flattened1^2)) * sqrt(sum(flattened2^2))
similarity <- dot_product / norm_product


sim300 <- similarity



results_simiralities <- compute_matrices_similarity(matrix_list)

length(iter_sets[["300"]])







# Check directories inside the data folder
habitats <- dir("data/performance_fastspar_iterations")
data_object <- list()

# Iterate over habitats
for (habitat in habitats) {
  habitat_path <- paste0("data/performance_fastspar_iterations/", habitat)

  # Get interactions within the habitat
  interactions <- dir(habitat_path)
  habitat_data <- list()

  # Iterate over interactions within the habitat
  for (interaction in interactions) {
    interaction_path <- file.path(habitat_path, interaction)
    tables <- dir(interaction_path, pattern = ".cor")
    interaction_data <- list()

    # Iterate over tables within the interaction
    for (table in tables) {
      table_path <- file.path(interaction_path, table)
      
      # Load the table data (adjust this step based on your specific requirements)
      table_data <- as.matrix(read_tsv(table_path)) #try remove as matrix
      
      #Store the table data within the interaction data
      interaction_data[[table]] <- table_data
     }
    
    # Store the interaction data within the habitat data
    habitat_data[[interaction]] <- interaction_data
  }
  
  # Store the habitat data within the data object
  data_object[[habitat]] <- habitat_data
}

