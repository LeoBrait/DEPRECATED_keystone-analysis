library("tidyverse")
source("R/data_processing/calculate_cosine_similarity.R")

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
      table_data <- read_tsv(table_path)
      
      # Store the table data within the interaction data
      #interaction_data[[table]] <- table_data
    }}}
    
#     # Store the interaction data within the habitat data
#     habitat_data[[interaction]] <- interaction_data
#   }
  
#   # Store the habitat data within the data object
#   data_object[[habitat]] <- habitat_data
# }

