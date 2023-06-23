options(scipen = 9999999)
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
matrix_list <- lapply(iter_sets$"16000", function(df) as.matrix(convert_to_numeric(df)))

x <- as.matrix(matrix_list[[1]])


x <- cosine_similarity(matrix_list[[1]], matrix_list[[2]])


flattened1 <- c(matrix_list[[1]])
flattened2 <- c(matrix_list[[2]])
dot_product <- sum(flattened1 * flattened2)
norm_product <- sqrt(sum(flattened1^2)) * sqrt(sum(flattened2^2))
similarity <- dot_product / norm_product


sim300 <- similarity



results_simiralities <- compute_matrices_similarity(matrix_list)

# get half of the matrix and convert to vector
results_simiralities <- as.vector(results_simiralities[upper.tri(results_simiralities, diag = FALSE)])


data <- as.data.frame(results_simiralities)

ggplot(data = data
        , aes(x = results_simiralities)) +
    geom_histogram(bins = 100) +
    labs(x = "Cosine similarity", y = "Frequency") +
    theme_bw(
)
ggsave("cosine_similarity.png", width = 10, height = 10, units = "cm")
