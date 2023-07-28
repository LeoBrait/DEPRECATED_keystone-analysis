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

args <- commandArgs(trailingOnly = TRUE)
analysis_frame <- as.character(args[1])
annotated_table_relative <- as.character(args[2])
metadata_table <- as.character(args[3])
parallel <- as.numeric(args[4])
minimum_samples <- as.numeric(args[5])

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


keystones <- read_csv(paste0(
    "data/", analysis_frame, "/final_keystones_table/keystones.csv")) %>%
  rename(habitat = Habitat) %>%
  mutate(habitat = gsub(" phyla", "", habitat)) %>%
  mutate(Taxon = gsub("\\.", " ", Taxon))

############################### Data Treatment #################################

phyla_abundances <- phyla_abundances %>%
  group_by(habitat) %>%
  mutate(n_samples = n()) %>%
  filter(n_samples >= minimum_samples) %>%
  ungroup() %>%
  select(-n_samples)

phyla_abundances <- phyla_abundances %>%
  filter(habitat %in% keystones$habitat)

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
      print("Output already exists!")
      load(paste0(results_path, "mdsgeral.RData"))
}

if (!file.exists(paste0(results_path, "permanova_ecosystem.RData"))) {
  print("Ecosystem's Permanova not found!, running...")
  permanova_ecosystem <- adonis2(
  community_distancematrix ~ ecosystem,
  data = phyla_abundances[1:8],
  permutations = 4999)
  save(
    permanova_ecosystem,
    file = paste0(results_path, "permanova_ecosystem.RData"))
    print("Done!")
    } else {
      print("Output already exists!")
      load(paste0(results_path, "permanova_ecosystem.RData"))
}

if (!file.exists(paste0(results_path, "permanova_lifestyle.RData"))) {
  print("Life-style's permanova not found!, running...")
  permanova_lifestyle <- adonis2(
  community_distancematrix ~ life_style,
  data = phyla_abundances[1:8],
  permutations = 4999)
  save(
    permanova_lifestyle,
    file = paste0(results_path, "permanova_lifestyle.RData"))
    print("Done!")
    } else {
      print("Output already exists!")
      load(paste0(results_path, "permanova_lifestyle.RData"))
}

if (!file.exists(paste0(results_path, "simper.RData"))) {
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
    print("Output already exists!")
    load(paste0(results_path, "simper.RData"))
}
