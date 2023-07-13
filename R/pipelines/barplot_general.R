source("R/src/draw_barplot_simple.R")
##1 LIFESTYLE ##################################################################
##1.1 RICHNESS
#Bonafide
bonafide_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "Bonafide"),
  title = "Bonafide",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 25, 50),
  break_labels = c("    0", "25", "50"),
  colors = life_style_colors)

#CPR
cpr_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "CPR"),
  title = "CPR",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 33, 65),
  break_labels = c("    0", "33", "65"),
  colors = life_style_colors)

#DPANN
dpann_barplot_richness_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "DPANN"),
  title = "DPANN",
  x_var = "life_style",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Life Style",
  legend_position = "none",
  breaks = c(0, 1.5, 3),
  break_labels = c("    0", "1.5", "3"),
  colors = life_style_colors)

## 1.2 ABUNDANCE
#bonafide
bonafide_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "Bonafide"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance (%)",
  title_x = "",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0, 0.50, 1),
  break_labels = c("    0", "50", "100"),
  colors = life_style_colors)

#CPR
cpr_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "CPR"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "Life style",
  legend_title = "Life style",
  legend_position = "none",
  breaks = c(0, 0.0060, 0.012),
  break_labels = c("    0", "0.6", "1.2"),
  colors = life_style_colors)

#DPANN
dpann_abundance_lifestyle <- draw_barplot_simple(
  data = subset(lifestyle_microgroups_prevalence, microGroup == "DPANN"),
  title = "",
  x_var = "life_style",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "Life style",
  legend_position = "top",
  breaks = c(0, 0.00012, 0.00023),
  break_labels = c("    0", "0.012", "0.023"),
  colors = life_style_colors)

## 2. ECOSSYSTEM ###############################################################
## 2.1 RICHNESS
#bonafide
bonafide_richness_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "Richness",
  title_x = "",
  legend_title = "Ecosystem",
  legend_position = "top",
  breaks = c(0, 29, 58),
  break_labels = c("    0", "29", "58"),
  colors = ecosystem_colors)

ecosystem_legend <- cowplot::get_legend(bonafide_richness_environment)
bonafide_richness_environment <- bonafide_richness_environment + 
  theme(legend.position = "none")

#cpr
cpr_richness_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Environment",
  legend_position = "none",
  breaks = c(0, 41, 83),
  break_labels = c("    0", "41", "83"),
  colors = ecosystem_colors)

#dpann
dpann_richness_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "DPANN"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_rich_life",
  error_var = "se_rich_life",
  title_y = "",
  title_x = "",
  legend_title = "Environment",
  breaks = c(0, 2.3, 4.6),
  break_labels = c("    0", "2.3", "4.6"),
  colors = ecosystem_colors)

## 2.2 ABUNDANCE
#bonafide
bonafide_abundance_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "Bonafide"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "Relative abundance (%)",
  title_x = "",
  legend_title = "",
  legend_position = "none",
  breaks = c(0, 0.5, 1),
  break_labels = c("    0", "50", "100"),
  colors = ecosystem_colors)

#cpr
cpr_abundance_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "CPR"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "Ecosystem",
  legend_title = "Environment",
  legend_position = "none",
  breaks = c(0, 0.06, 0.12),
  break_labels = c("    0", "6", "12"),
  colors = ecosystem_colors)

#dpann
dpann_abundance_environment <- draw_barplot_simple(
  data = subset(ecosystem_microgroups_prevalence, microGroup == "DPANN"),
  title = "",
  x_var = "ecosystem",
  y_var = "mean_abu_life",
  error_var = "se_abu_life",
  title_y = "",
  title_x = "",
  legend_title = "Environment",
  legend_position = "none",
  breaks = c(0, 0.0010, 0.00215),
  break_labels = c("    0", "0.10", "0.21"),
  colors = ecosystem_colors)


### GET LEGENDS ################################################################

life_style_legend <- cowplot::get_legend(dpann_abundance_lifestyle)
 dpann_abundance_lifestyle <- dpann_abundance_lifestyle +
  theme(legend.position = "none")

######### MERGE PLOTS ##########################################################
barplot_abundance_lifestyle <- ggarrange(
  bonafide_abundance_lifestyle,
  cpr_abundance_lifestyle,
  dpann_abundance_lifestyle,
  ncol = 3)

barplot_richness_lifestyle <- ggarrange(
  bonafide_barplot_richness_lifestyle,
  cpr_barplot_richness_lifestyle,
  dpann_barplot_richness_lifestyle,
  ncol = 3)

barplot_richness_environment <- ggarrange(
  bonafide_richness_environment,
  cpr_richness_environment,
  dpann_richness_environment,
  ncol = 3)

barplot_abundance_environment <- ggarrange(
  bonafide_abundance_environment,
  cpr_abundance_environment,
  dpann_abundance_environment,
  ncol = 3)