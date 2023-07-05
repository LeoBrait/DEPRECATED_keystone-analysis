library("tidyverse")
library("stringr")
library("cowplot")
source("R/src/draw_lineplot.R")

args <- commandArgs(trailingOnly = TRUE)
frame_analysis <- as.character(args[1])


if (!file.exists("results/")) {
  dir.create("results/")
  }
if (!file.exists(paste0("results/", frame_analysis))) {
  dir.create(paste0("results/", frame_analysis))
  }

summary_path <- paste0(
    "data/", frame_analysis, "/summaries/performance_fastspar_iterations/")

plot_path <- paste0(
    "results/", frame_analysis, "/performance_fastspar_iterations/")

# orderer
animal_feces_df <- read_csv(paste0(
    summary_path, "animal_host-associated.animal_feces.csv")) %>%
    gather("iteration", "cosine_similarity")
orderer <- as.numeric(str_extract(animal_feces_df$iteration, "\\d+"))
orderer <- unique(orderer)
orderer <- as.factor(orderer)


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

groundwatermine_df <- read_csv(paste0(
    summary_path, "groundwater.mine.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

groundwaterporous_df <- read_csv(paste0(
    summary_path, "groundwater.porous_contaminated.csv")) %>%
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

soiltundra_df <- read_csv(paste0(
    summary_path, "soil.tundra_soil.csv")) %>%
    gather("iteration", "cosine_similarity") %>%
    mutate(iteration = as.factor(iteration)) %>%
    mutate(iteration = factor(iteration, levels = orderer))

soiltundra_plot <- draw_lineplot(
    data = soiltundra_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Soil tundra (N = 3)",
    x_title = " ",
    y_title = "Cosine similarity (%)") +
    scale_x_discrete(labels = NULL) +
    geom_vline(xintercept = "4000", linetype = "dashed")

groundwatermine_plot <- draw_lineplot(
    data = groundwatermine_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Groundwater mine (N = 3)",
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

groundwaterporous_plot <- draw_lineplot(
    data = groundwaterporous_df,
    x_var = "iteration",
    y_var = "cosine_similarity",
    habitat_name = "Groundwater porous (N = 48)",
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
    soiltundra_plot,
    groundwatermine_plot,
    aqueoushumour_plot,
    salinehyper_plot,
    soilsavanna_plot,
    groundwaterporous_plot,
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

