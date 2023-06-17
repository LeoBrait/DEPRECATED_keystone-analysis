
source("R/src/utilities/install_and_load.R")
install_and_load(libs = c("tidyverse" = "2.0.0"))



source("R/src/data_wrangling/merge_annotation_metadata.R")
taxon_abundances <- merge_annotation_metadata(
    annotation_df = read_csv(
        "data/input/taxon_abundances/kraken_custom_phyla.csv"),
    metadata_df = read_csv("data/input/biome_classification.csv"),
    metadata_variables = c("samples", "habitat")
    )

habitat_up12 <- taxon_abundances %>%
    group_by(habitat) %>%
    summarise(n = n()) %>%
    filter(n > 12) %>%
    pull(habitat)


taxon_abundances <- taxon_abundances %>%
    filter(habitat %in% habitat_up12) %>%
    mutate(habitat = factor(habitat))


for (i in 1:length(habitat_up12)){
  subset <- subset(taxon_abundances, habitat == habitat_up12[i])
  subset <- subset %>% select(-habitat)

  #getting rid of zeros
  subset_numeric <- subset[, 2:ncol(subset)]
  subset_numeric_clean <- subset_numeric[, colSums(subset_numeric != 0) > 0]
  subset_clean <- cbind(samples = subset$samples, subset_numeric_clean)

  write.csv(
    subset_clean,
    paste0("community_matrix/", habitat_up12[i], ".phyla", ".csv"), row.names = FALSE)
}
