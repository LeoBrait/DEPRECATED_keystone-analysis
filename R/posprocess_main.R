# title: Keystones Paper
# authors:
#  Camilo M. Ferreira,
#  Leonardo Brait,
#  Pablo Viana,
#  Felipe Alexandre
################################## Environment #################################
set.seed(140822)
options(scipen = 9999999)

if (interactive()) {
  analysis_frame <- "phyla_analysis_aug21"
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

source("R/src/install_and_load.R")
install_and_load(
  libs = c(
    "logging" = "0.10-108",
    "janitor" = "2.2.0",
    "lmPerm" = "2.1.0",
    "multcomp" = "1.4-24",
    "sf" = "1.0-13",
    "cowplot" = "1.1.1",
    "maps" = "3.4.1",
    "ggmap" = "3.0.2",
    "ggpubr" = "0.6.0",
    "patchwork" = "1.1.2",
    "rnaturalearth" = "0.3.3",
    "rnaturalearthdata" = "0.1.0",
    "ggExtra" = "0.10.0",
    "svglite" = "2.1.1",
    "mgcv" = "1.8-42",
    "mgcViz" = "0.1.9",
    "gdata" = "2.19.0",
    "gridExtra" = "2.3",
    "emmeans" = "1.8.6",
    "vegan" = "2.6-4",
    "ggpubr" = "0.6.0",
    "ggh4x" = "0.2.4",
    "funrar" = "1.5.0",
    "multcompView" = "0.1-9",
    "rpart" = "4.1.19",
    "ggh4x" = "0.2.4",
    "reshape2" = "1.4.4",
    "scales" = "1.2.1",
    "randomForest" = "4.7-1.1",
    "tidyverse" = "2.0.0",
    "hrbrthemes" = "any"
  )
)

################################## Load data ###################################
source("R/src/merge_annotation_metadata.R")

phyla_abundances <-
  merge_annotation_metadata(
    annotation_df = read.csv(
      paste0("data/taxon_abundances/", annotated_table_relative)
    ),
    metadata_df = read.csv(
      paste0("data/metadata/", metadata_table)
    ),
    metadata_variables = c(
      "samples",  "biosphere", "ecosystem",
      "habitat", "life_style",  "latitude",
      "longitude"
    )
  )

network_scores <-
  read_csv(
    paste0("data/", analysis_frame, "/final_keystones_table/keystones.csv")
  ) %>%
  rename(habitat = Habitat) %>%
  rename(ecosystem = Ecosystem)

radiation <- read_csv("data/radiations/radiations_phyla.csv")

source("R/pipelines/checking_data_integrity.R")
############################### Data Treatment #################################

phyla_abundances <- phyla_abundances %>%
  group_by(ecosystem, habitat) %>%
  mutate(n_samples = n()) %>%
  filter(n_samples >= 5) %>%
  ungroup() %>%
  select(-n_samples)


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

############################### Visual settings ################################
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

########################## Validating custom database ##########################

#source("R/pipelines/validating_databases.R")

################## Data Wrangling (DPANN, CPR, BONAFIDE) #######################
# Description: Identify DPAAN, CPR and Bonafide in the samples and sites. Get
# the abundance and richness of these.

source("R/pipelines/obtaining_richness_abundance.R")

################################# Load Rdata ###################################

results_path <- paste0("data/", analysis_frame, "/rdata/")
load(paste0(results_path, "mdsgeral.RData"))
load(paste0(results_path, "permanova_ecosystem.RData"))
load(paste0(results_path, "permanova_lifestyle.RData"))

################################# Plot Pannels #################################

source("R/pipelines/drawing_panel1.R")
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