# title: Keystones Paper
# authors:
#  Camilo M. Ferreira,
#  Leonardo Brait,
#  Pablo Viana,
#  Felipe Alexandre
################################## Environment #################################

if (!file.exists("r_libs")) {
  dir.create("r_libs")
}

source("R/src/install_and_load.R")
install_and_load(
  libs = c(
    "tidyverse" = "1.2.1",
    "vegan" = "2.5-2"
  ),
  loc = "r_libs"
)

if (interactive()) {
  analysis_frame <- "phyla_analysis_july23"
  annotated_table_relative <- "annotated_table_relative.csv"
  annotated_table_relative <- "kraken_biomedb_relative_phyla.csv"
  metadata_table <- "biome_classification.csv"
  parallel <- 40
  minimum_samples <- 5
} else {
  args <- commandArgs(trailingOnly = TRUE)
  analysis_frame <- as.character(args[1])
  annotated_table_relative <- as.character(args[2])
  metadata_table <- as.character(args[3])
  parallel <- as.numeric(args[4])
  minimum_samples <- as.numeric(args[5])
}

set.seed(140822)
options(scipen = 9999999)

################################## Load data ###################################
source("R/src/merge_annotation_metadata.R")
phyla_abundances <- merge_annotation_metadata(
    annotation_df = read.csv(paste0(
        "data/taxon_abundances/", annotated_table_relative)),
    metadata_df = read.csv(paste0(
        "data/metadata/", metadata_table)),
    metadata_variables = c(
        "samples",
        "biosphere",
        "ecosystem",
        "habitat",
        "life_style",
        "latitude",
        "longitude"))

############################### Data Treatment #################################

phyla_abundances <- phyla_abundances %>%
  group_by(ecosystem, habitat) %>%
  mutate(n_samples = n()) %>%
  filter(n_samples >= minimum_samples) %>%
  ungroup() %>%
  select(-n_samples)

# Visual treatment
source("R/src/visual_treat.R")
phyla_abundances <- treatment(phyla_abundances)

############## nmds analysis and permanova for entire community ################

results_path <- paste0("data/", analysis_frame, "/rdata/")
if (!file.exists(results_path)) {
  dir.create(results_path)
}

standarized_phyla <- decostand(phyla_abundances[8:length(phyla_abundances)],
method = "hellinger")

community_distancematrix <- vegdist(
  x = standarized_phyla,
  method = "jaccard")


if (!file.exists(paste0(results_path, "mdsgeral.RData"))) {
  print("MDS not found!, running... This can take a lot of time!")
  ord <- metaMDS(
  standarized_phyla,
  distance = "jaccard",
  try = 1, trymax = 4999, stress = 1, parallel = parallel)
  save(
    ord,
    file = paste0(results_path, "mdsgeral.RData"))
    print("Done!")
    } else {
      print("MDS Output already exists!")
      load(paste0(results_path, "mdsgeral.RData"))
}

if (!file.exists(paste0(results_path, "permanova_ecosystem.RData"))) {
  print("Ecosystem's Permanova not found!, running...")
  permanova_ecosystem <- adonis2(
  community_distancematrix ~ ecosystem,
  data = phyla_abundances[1:8],
  permutations = 4999,
  parallel = parallel)
  save(
    permanova_ecosystem,
    file = paste0(results_path, "permanova_ecosystem.RData"))
    print("Done!")
    } else {
      print("Ecosystem's Permanova output already exists!")
      load(paste0(results_path, "permanova_ecosystem.RData"))
}

if (!file.exists(paste0(results_path, "permanova_lifestyle.RData"))) {
  print("Life-style's permanova not found!, running...")
  permanova_lifestyle <- adonis2(
  community_distancematrix ~ life_style,
  data = phyla_abundances[1:8],
  permutations = 4999,
  parallel = parallel)
  save(
    permanova_lifestyle,
    file = paste0(results_path, "permanova_lifestyle.RData"))
    print("Done!")
    } else {
      print("Life-style's permanova output already exists!")
      load(paste0(results_path, "permanova_lifestyle.RData"))
}


# cluster <- phyla_abundances %>%
#   group_by(ecosystem, habitat) %>%
#   mutate(cluster = paste0(ecosystem, "_", habitat)) %>%
#   pull(cluster)

if (!file.exists(paste0(results_path, "simper.RData"))) {
  print("Simper not found!, running...")
  category <- phyla_abundances[["ecosystem"]]
  raw_simper <-  simper(
    standarized_phyla,
    group = category,
    parallel = parallel,
    permutations = 4999
  )
  save.image(paste0(results_path, "simper.RData"))
  print("Done!")
  } else {
    print("Simper output already exists!")
    load(paste0(results_path, "simper.RData"))
}

summary_simper <- summary(raw_simper)

for (name_of_table in names(summary_simper)){
  df <- summary_simper[[name_of_table]]
  prev <- 0

  for (rows in 1:nrow(df)){
    df[rows, 8] <-  df[rows, 6] - prev
    prev <-  df[rows, 6]
  }

  df$comparison <- rep(name_of_table, nrow(df))
  colnames(df) <- c(
      "average",            "sd",     "ratio",
          "ava",           "avb",    "cumsum",
            "p",  "contribution", "comparison"
  )
  df <- tibble::rownames_to_column(df, "OTU")
  simpertables[[name_of_table]] <- df
}

#Reuniting all data in a single data.frame
simper_clean <- bind_rows(simpertables)
simper_clean <- simper_clean[simper_clean$p < 0.05, ]

write.csv(
  file = paste0(results_path, "summaries/simper_ecosystem.csv"),
  x = simper_clean
)
