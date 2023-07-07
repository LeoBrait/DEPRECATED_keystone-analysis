# title: Keystones Paper
# authors:
#  Camilo M. Ferreira,
#  Leonardo Brait,
#  Pablo Viana,
#  Felipe Alexandre
################################## Environment #################################
library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)
analysis_frame <- as.character(args[1])
annotated_table_relative <- as.character(args[2])
metadata_table <- as.character(args[3])
parallel <- as.numeric(args[4])

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


radiation <- read_csv(paste0(
    "data/radiations/", analysis_frame, "/radiations.csv"))

source("R/pipelines/checking_data_integrity.R")


############################### Data Treatment #################################

phyla_abundances <- phyla_abundances %>%
  filter(habitat %in% keystones$habitat)

dpann_groups <- radiation %>%
  filter(microgroup == "DPANN") %>%
  pull(taxon)

cpr_groups <- radiation %>%
  filter(microgroup == "CPR") %>%
  pull(taxon)

phyla_abundances_long <- phyla_abundances %>%
  gather(
    Taxon,
    abundance,
    -samples,
    -biosphere,
    -ecosystem,
    -habitat,
    -life_style,
    -latitude,
    -longitude) %>%
  mutate(Taxon = gsub("\\.", " ", Taxon)) #TODO: remove

phyla_abundances_habitat <- phyla_abundances_long %>%
  group_by(life_style, ecosystem, habitat, Taxon) %>%
  summarise(
        mean_abundance = mean(abundance),
        sd_abundance = sd(abundance),
        se_abundance = sd_abundance / sqrt(n())) %>%
  ungroup()


keystones_abundance_habitat <- phyla_abundances_habitat %>%
 left_join(keystones, by = c("Taxon", "habitat")) %>%
 mutate(
    microGroup = factor(case_when(
      Taxon %in% dpann_groups ~ "DPANN",
      Taxon %in% cpr_groups ~ "CPR",
      TRUE ~ "Bonafide"))) %>%
      complete(Taxon) %>%
  replace_na(list(
    EDpDM_isKeystone = 0,
    AbundanceAbsolute_zScore = 0,
    EDpDM_zScore = 0,
    BC_zScore = 0,
    LIASPdir_zScore = 0,
    LIASPindir_zScore = 0,
    Abundance = 0,
    LIASP_isKeystone = 0,
    LIASP = 0,
    LIASPDir_zScore = 0,
    LIASPindir  = 0,
    EDpDM_rank = 0,
    mean_abundance = 0))


# Visual treatment
source("R/src/visual_treat.R")
phyla_abundances <- treatment(phyla_abundances)
phyla_abundances_long <- treatment(phyla_abundances_long)
phyla_abundances_habitat <- treatment(phyla_abundances_habitat)
keystones_abundance_habitat <- treatment(keystones_abundance_habitat)

# ############################## Visual settings #################################
ecosystem_colors <- c(
    "Human Host" = "#DA70D6",
    "Animal Host" = "#FFD700",
    "Plant Host" = "#228B22",
    "Groundwater" = "#3A5FCD",
    "Freshwater" = "#87CEFA",
    "Wastewater" = "#000000",
    "Saline Water" = "#20B2AA",
    "Sediment" = "#F4A460",
    "Soil" = "#8B4513")

life_style_colors <- c(
    "Holobiont" = "#ED1A1A",
    "Free-living" = "#012169")

microgroup_colors <- c(
  "Bonafide" = "#009739",
  "CPR" = "#FEDD00",
  "DPANN" = "#012169"
                        )

# ######################## Validating custom database ##########################

# source("source//pipelines//validating_databases.R")

# ############## Data Wrangling (DPANN, CPR, BONAFIDE) #########################
# # Description: Identify DPAAN, CPR and Bonafide in the samples and sites. Get
# # the abundance and richness of these.

# source("source/pipelines/obtaining_richness_abundance.R")

############## nmds analysis and permanova for entire community ################
# WARNING: please extract the designated zip files in the outputs directory
# before run this code. It is very expensive.
results_path <- paste0("results/", analysis_frame, "/rdata/")
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
    } else {
      print("Output already exists!")
      load(paste0(results_path, "mdsgeral.RData"))
  }

if (!file.exists(paste0(results_path, "permanova_ecosystem.RData"))) {
  print("Ecosystem's Permanova not found!, running...")
  permanova <- adonis2(
  community_distancematrix ~ ecosystem + habitat + life_style,
  data = phyla_abundances[1:8],
  permutations = 4999)
  save(
    permanova,
    file = paste0(results_path, "permanova_ecosystem.RData"))
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
    } else {
      print("Output already exists!")
      load(paste0(results_path, "permanova_lifestyle.RData"))
  }

category <- phyla_abundances[["ecosystem"]]
raw_simper <-  simper(
                  standarized_phyla,
                  group = category,
                  parallel = parallel,
                  permutations = 4999
                      )
save.image(paste0(results_path, "simper.RData"))


# ################################# Plot Pannels #################################

# source("source//panels//panel1.R")
# source("source//panels//panel2.R")
# source("source//panels//panel3.R")
# source("source//panels//panel4.R")

# ############################ General Additive Models ###########################
# source("source/analysis/do_gam.R")
# bonafide_richnessmodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "Bonafide",
#   response = "meanRich",
#   k = 20,
#   output = "bonafide_richness")
# bonafide_abundancemodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "Bonafide",
#   response = "meanAbu",
#   k = 20,
#   output = "bonafide_abundance")

# cpr_richnessmodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "CPR",
#   response = "meanRich",
#   k = 20,
#   output = "cpr_richness")
# cpr_abundancemodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "CPR",
#   response = "meanAbu",
#   k = 20,
#   output = "cpr_abundance")

# dpann_richnessmodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "DPANN",
#   response = "meanRich",
#   k = 20,
#   output = "dpann_richness")
# dpann_abundancemodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = "DPANN",
#   response = "meanAbu",
#   k = 20,
#   output = "dpann_abundance")

# complete_richnessmodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = c("Bonafide", "CPR", "DPANN"),
#   response = "meanRich",
#   k = 20,
#   output = "complete_richness")
# complete_abundancemodel <- do_gam(data = microgroups_prevalence_persite,
#   microgroup = c("Bonafide", "CPR", "DPANN"),
#   response = "meanAbu",
#   k = 20,
#   output = "complete_abundance")