#' @title Ploting iterations performance
#' @description This script plots the performance of the iterations of the
#' fastspar algorithm based on the cosine similarity between the matrices.
#' @author Brait LAS.
#' @date july 2023.
#' @reviewer None.

################################# Environment ##################################
options(scipen = 9999999)
options(digits = 20)

source("R/src/install_and_load.R")
install_and_load(
  libs = c(
    "stringr" = "any",
    "ggpubr" = "0.2.4",
    "cowplot" = "0.9.4",
    "tidyverse" = "any"
  )
)

if (interactive()) {
  analysis_frame <- "phyla_analysis_july23"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  analysis_frame <- as.character(args[1])
}

if (!file.exists("results/")) {
  dir.create("results/")
}
if (!file.exists(paste0("results/", analysis_frame))) {
  dir.create(paste0("results/", analysis_frame))
}

summary_path <-
    paste0(
        "data/", analysis_frame, "/summaries/performance_fastspar_iterations/"
    )

plot_path <-
    paste0(
        "results/", analysis_frame, "/performance_fastspar_iterations/"
    )

####################### Reading and preparing the data #########################

# ordering the iterations --------------
animal_feces_df <- read_csv(
    paste0(summary_path, "animal_host-associated.animal_feces.csv")) %>%
    gather("iteration", "cosine_similarity")
orderer <- as.numeric(str_extract(animal_feces_df$iteration, "\\d+"))
orderer <- unique(orderer)
orderer <- as.factor(orderer)
orderer <- levels(orderer)


# reading the data ---------------------
animal_feces_df <- read_csv(paste0(
    summary_path, "animal_host-associated.animal_feces.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

aqueous_humour_df <- read_csv(paste0(
    summary_path, "animal_host-associated.aqueous_humour.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

groundwater_porouscont_df <- read_csv(paste0(
    summary_path, "groundwater.porous_contaminated.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

sulfurspring_df <- read_csv(paste0(
    summary_path, "freshwater.sulfur_spring.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

humangut_df <- read_csv(paste0(
    summary_path, "human_host-associated.human-gut.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

salinecoastal_df <- read_csv(paste0(
    summary_path, "saline_water.coastal_seawater.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

saline_hyper_df <- read_csv(paste0(
    summary_path, "saline_water.hypersaline_water.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

soilsavanna_df <- read_csv(paste0(
    summary_path, "soil.savanna_soil.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

sediment_lakesed_df <- read_csv(paste0(
    summary_path, "sediment.lake_sediment.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

################################## Plotting ####################################
source("R/src/draw_lineplot.R")


lakesed_plot <- draw_lineplot(
    data = sediment_lakesed_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Lake Sediment (N = 5)",
    x_title = " ",
    y_title = "Cosine similarity (%)") +
    scale_x_discrete(labels = NULL) +
    geom_vline(xintercept = "4000", linetype = "dashed")

sulfurspring_plot <- draw_lineplot(
    data = sulfurspring_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Sulfur Spring (N = 6)",
    x_title = " ",
    y_title = " ") +
    scale_x_discrete(labels = NULL) +
    geom_vline(xintercept = "5000", linetype = "dashed")

aqueoushumour_plot <- draw_lineplot(
    data = aqueous_humour_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Aqueous humour (N = 8)",
    x_title = " ",
    y_title = " ") +
    scale_x_discrete(labels = NULL)  +
    geom_vline(xintercept = "4000", linetype = "dashed")

salinehyper_plot <- draw_lineplot(
    data = saline_hyper_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Saline hypersaline (N = 16)",
    x_title = " ",
    y_title = "Cosine similarity (%)") +
    scale_x_discrete(labels = NULL)  +
    geom_vline(xintercept = "4000", linetype = "dashed")

soilsavanna_plot <- draw_lineplot(
    data = soilsavanna_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Soil savanna (N = 21)",
    x_title = " ",
    y_title = " ") +
    scale_x_discrete(labels = NULL)  +
    geom_vline(xintercept = "4000", linetype = "dashed")

groundwater_porouscont_plot <- draw_lineplot(
    data = groundwater_porouscont_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Groundwater porous Contaminated (N = 48)",
    x_title = " ",
    y_title = " ")  +
    scale_x_discrete(labels = NULL)  +
    geom_vline(xintercept = "4000", linetype = "dashed")

humangut_plot <- draw_lineplot(
    data = humangut_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Human gut (N = 58)",
    x_title = "Iterations",
    y_title = "Cosine similarity (%)")  +
    geom_vline(xintercept = "4000", linetype = "dashed")

salinecoastal_plot <- draw_lineplot(
    data = salinecoastal_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Saline coastal (N = 286)",
    x_title = "Iterations",
    y_title = " ") +
    geom_vline(xintercept = "4000", linetype = "dashed")

animalfeces_plot <- draw_lineplot(
    data = animal_feces_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Animal feces (N = 675)",
    x_title = "Iterations",
    y_title = " ") +
    geom_vline(xintercept = "4000", linetype = "dashed")

pannel <- plot_grid(
    lakesed_plot,
    sulfurspring_plot,
    aqueoushumour_plot,
    salinehyper_plot,
    soilsavanna_plot,
    groundwater_porouscont_plot,
    humangut_plot,
    salinecoastal_plot,
    animalfeces_plot,
    ncol = 3)

ggsave(
    filename = paste0(
        plot_path, "cosine_similarity_lineplot.png"),
    plot = pannel,
    width = 12,
    height = 8,
    dpi = 300)
ggsave(
    filename = paste0(
        plot_path, "cosine_similarity_lineplot.svg"),
    plot = pannel,
    width = 12,
    height = 8,
    dpi = 300)
