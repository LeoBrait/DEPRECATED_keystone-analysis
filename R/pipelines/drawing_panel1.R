#####  World Map of samples  ###################################################
library("rnaturalearth")
library("sf")
grouped_samples <- phyla_abundances_long %>%
  group_by(latitude, longitude, ecosystem) %>%
  summarise(nSamples = n()) %>%
  ungroup() %>%
  mutate(
    nSamples = as.numeric(nSamples),
    gradient = case_when(
      latitude < -66.5 ~ "South Pole",
      latitude > 66.5 ~ "North Pole",
      latitude >= -66.5 & latitude <= -35.5 ~ "South Temperate",
      latitude >= -35.5 & latitude <= -23.5 ~ "South Subtropical",
      latitude >= -23.5 & latitude <= 0 ~ "South Tropical",
      latitude >= 0 & latitude <= 23.5 ~ "North Tropical",
      latitude >= 23.5 & latitude <= 35.5 ~ "North Subtropical",
      latitude >= 35.5 & latitude <= 66.5 ~ "NorthTemperate"))

worldmap <- ggplot(ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf(fill = "#f9f7f0") +
  coord_sf(expand = FALSE) +

  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +

  geom_point(data = grouped_samples,
    aes(
      x = longitude,
      y = latitude,
      size = nSamples,
      color = ecosystem),
    alpha = 0.6,
    shape = 19) +
  scale_color_manual(values = ecosystem_colors)

########## Richness vs Latitude ################################################
source("R/src/draw_latitude_gam.R")
bonafide_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microGroup == "Bonafide"),
  breaks = c(45, 55, 65),
  title_y = "Richness",
  main_title = "Bonafide",
  title_x = "")

cpr_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microGroup == "CPR"),
  breaks = c(45, 68, 84.5),
  title_y = "",
  main_title = "CPR",
  title_x = "Latitude")

dpann_latitude_plot <- draw_latitude_gam(
  data = microgroups_prevalence_persite %>% filter(microGroup == "DPANN"),
  breaks = c(2.10, 3.55, 4.3),
  title_y = "",
  main_title = "DPANN",
  title_x = "")

plot_latitude <- plot_grid(
  bonafide_latitude_plot,
  cpr_latitude_plot,
  dpann_latitude_plot,
  ncol = 3,
  rel_widths = c(1, 1, 1))

############# barplots for richness and abundance ##############################
source("R/pipelines/barplot_general.R")

############# nmds plot #######################################################

ord_df <- cbind(
  phyla_abundances$ecosystem,
  as.data.frame(ord$points),
  phyla_abundances$samples) %>%
  mutate(stress = ord$stress)
colnames(ord_df) <- c("ecosystem", "MDS1", "MDS2", "samples", "stress")

##Plotting nmds
nmds_allsamples <- ggplot(
  ord_df,
  aes(x = MDS1, y = MDS2, color = ecosystem)) +

  ## Theme
  theme_pubr() +
  theme(
    text = element_text(size = unit(9, "cm"), family = "sans"),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.margin = margin(t = 3.5, unit = "cm"),
    legend.text = element_text(size = 9, family = "sans"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(size = unit(9, "cm"), family = "sans")#,
  ) +

  ## Points
  geom_point(size = 1.2, shape = 20) +
  scale_color_manual(values = ecosystem_colors, name = "Ecosystem") +

  ##axis text
  theme(axis.text = element_text(face = "bold", size = unit(9, "cm"))) +
  ggtitle("") +
  labs(
    x = paste0(
     "MDS1 (", round(attr(ord$species, "shrinkage")[1] * 100,
     digits = 2), " %)"
              ),
    y = paste0(
      "MDS2 (", round(attr(ord$species, "shrinkage")[2] * 100,
      digits = 2), " %)"
              )
      ) +

  ##Anotation
  annotate("text",
    x = min(ord_df$MDS1) + 2,
    y = max(ord_df$MDS2) - 0.5,
    label = paste(
      "Stress =",
      round(unique(ord_df$stress)[1],
        digits = 3),
      "\n",
      "Life style R² =", round(unique( permanova_lifestyle$R2)[1],
        digits = 3),
      "\n",
      "p-value < ", round(unique( permanova_lifestyle$"Pr(>F)")[1],
        digits = 4),
      "\n",
      "Ecossystem R² =", round(unique(permanova$R2)[1],
        digits = 3),
      "\n",
      "p-value <", round(unique(permanova$"Pr(>F)")[1],
        digits = 4)
    ),
    size = unit(3, "cm")) +

  ##Legend
    guides(color = guide_legend(override.aes = list(size = 4, shape = 16))) +
    annotation_custom(
      life_style_legend,
      xmin = 1,
      xmax = 3.1,
      ymin = 1,
      ymax = 2) +

  ## breaks
  scale_y_continuous(breaks = c(-1, 0, 1, 2)) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2)) +
  ylim(-1, 2)

####################### Mergering plots ########################################

top_right <- plot_grid(
  barplot_richness_lifestyle, barplot_abundance_lifestyle,
  ncol = 1, align = "hv",
  rel_widths = c(1, 1),
  labels = c("D", "E"))

bottom_right <- plot_grid(
  barplot_richness_environment, barplot_abundance_environment,
  ncol = 1, align = "hv",
  rel_heights = c(1, 1),
  labels = c("F", "G"))

right_panel <- plot_grid(
  top_right,
  bottom_right,
  nrow = 2,
  rel_widths = c(1, 1))

left_panel <- plot_grid(worldmap, nmds_allsamples, plot_latitude,
  ncol = 1,
  rel_heights = c(1.5, 2, 1),
  labels = c("A", "B", "C"))

panel1 <- plot_grid(left_panel, right_panel, ncol = 2, rel_widths = c(1, 1))

result_dir <- paste0("results/", analysis_frame, "/panel1/")
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

ggsave(
    paste0(result_dir, "panel1.pdf"),
    plot = panel1,
    width = 11,
    height = 8)

ggsave(
    paste0(result_dir, "panel1.png"),
    plot = panel1,
    width = 11,
    height = 8)

ggsave(
    paste0(result_dir, "panel1.svg"),
    plot = panel1,
    width = 11,
    height = 8)