# load standard database
phyla_abundances_stdb <- merge_annotation_metadata(
    annotation_df = read.csv(paste0(
        "data/taxon_abundances/", annotated_table_relative_stdb)),
    metadata_df = read.csv(paste0(
        "data/metadata/", metadata_table)),
    metadata_variables = c(
        "samples",
        "biosphere",
        "ecosystem",
        "habitat",
        "life_style",
        "latitude",
        "longitude")) %>%
    treatment()

phyla_abundances_long_std <- phyla_abundances_std %>%
  gather(
    Taxon,
    abundance,
    -samples,
    -biosphere,
    -ecosystem,
    -habitat,
    -life_style,
    -latitude,
    -longitude)

################################ RICHNESS ######################################
richness_std <- phyla_abundances_std %>%
    mutate(std_richness = rowSums(.[8:length(.)] > 0)) %>%
    select(samples, std_richness)

richness_custom <- phyla_abundances %>%
    mutate(custom_richness = rowSums(.[8:length(.)] > 0)) %>%
    select(samples, custom_richness) %>%
    filter(samples %in% richness_std$samples) # this shoudn't be necessary

phyla_richness <- inner_join(
        richness_std,
        richness_custom,
        by = c("samples")) %>%
    mutate(
        richness_diff = custom_richness - std_richness,
        richness_diff_percentage = richness_diff / std_richness * 100) %>%
    merge_annotation_metadata(
        metadata_df = read.csv(paste0("data/metadata/", metadata_table)),
        metadata_variables = c(
            "samples",
            "biosphere",
            "ecosystem",
            "habitat",
            "life_style",
            "latitude",
            "longitude")) %>%
    treatment()


#Summarize richness by ecosystem
phyla_richness_summ <- phyla_richness %>%
    group_by(ecosystem) %>%
    summarise(richness_diff_percentage = mean(richness_diff_percentage),
    richness_sd = sd(.$richness_diff_percentage),
    richness_se = richness_sd / sqrt(n()))

#################### SEQUENCES #################################################
sequences_custom <- read_csv("inputs/unclassified_count.csv")

sequences_std <- read_csv("inputs/standard_unclassified_count.csv")

sequences_comparison <- sequences_std %>%
    inner_join(sequences_custom,
        by = "samples",
        suffix = c("_stddb", "_customdb")) %>%
    merge_annotation_metadata(
        metadata_df = read.csv("inputs/final_biome_classification.csv"),
        metadata_variables = c(
            "samples",
            "biosphere",
            "ecosystem",
            "habitat",
            "life_style",
            "latitude",
            "longitude")) %>%
    treatment() %>%
    mutate(
            unclassified_diff =   unclassified_percentage_stddb
                                - unclassified_percentage_customdb)


sequences_comparison_summ <- sequences_comparison %>%
    group_by(ecosystem) %>%
    summarise(unclassified_diff = mean(unclassified_diff),
    unclassified_diff_sd = sd(.$unclassified_diff),
    unclassified_diff_se = unclassified_diff_sd / sqrt(n()))

#################### MICROGROUPS ###############################################

candidate_richness_std <- phyla_abundances_long_std %>%
    mutate(
        microGroup = factor(case_when(
        Taxon %in% dpann_groups ~ "DPANN",
        Taxon %in% cpr_groups ~ "CPR",
        TRUE ~ "Bonafide"
        ))) %>%
    group_by(
        samples,
        life_style,
        biosphere,
        ecosystem,
        habitat,
        latitude,
        longitude,
        microGroup) %>%
    filter(abundance != 0) %>%
    summarise(richness = n_distinct(Taxon)) %>%
    complete(microGroup) %>%
    replace_na(list(richness = 0)) %>%
    filter(!microGroup == "Bonafide") %>%
    summarise(richness = sum(richness))

candidate_richness_custom <- phyla_abundances_long %>%
  mutate(
    microGroup = factor(case_when(
      Taxon %in% dpann_groups ~ "DPANN",
      Taxon %in% cpr_groups ~ "CPR",
      TRUE ~ "Bonafide"
    ))) %>%
  group_by(
    samples,
    life_style,
    biosphere,
    ecosystem,
    habitat,
    latitude,
    longitude,
    microGroup) %>%
  filter(abundance != 0) %>%
  summarise(richness = n_distinct(Taxon)) %>%
  complete(microGroup) %>%
  replace_na(list(richness = 0)) %>%
  filter(!microGroup == "Bonafide") %>%
  summarise(richness = sum(richness)) %>%
  filter(samples %in% candidate_richness_std$samples) # this shoudn't be necessary #nolint

comparison_microgroups <- candidate_richness_custom %>%
    inner_join(
        candidate_richness_std,
        by = c("samples",
               "life_style",
               "biosphere",
               "ecosystem",
               "habitat",
               "latitude",
               "longitude"),
        suffix = c("_custom", "_std")) %>%
    mutate(
        richness_diff = richness_custom - richness_std,
        richness_diff_percentage = richness_diff / richness_std * 100)

#Summarize richness by ecosystem
richness_summary <- comparison_microgroups %>%
    group_by(ecosystem) %>%
    summarize(mean_richness = mean(richness_diff),
    sd_richness = sd(richness_diff),
    se_richness = sd_richness / sqrt(n()))

####### PLOTS ##################################################################
source("source/plots/draw_barplot_simple.R")

rich_gain_plot <- draw_barplot_simple(
    data = phyla_richness_summ,
    x_var = "ecosystem",
    y_var = "richness_diff_percentage",
    color = rep("gray", 9),
    error_var = "richness_se",
    title_y = "Gain in richness (%)",
    breaks = c(0, 50, 100, 150, 200, 250),
    break_labels = c("0", "50", "100", "150", "200", "250"))


classification_gain_plot <- draw_barplot_simple(
    data = sequences_comparison_summ,
    x_var = "ecosystem",
    y_var = "unclassified_diff",
    color = rep("gray", 9),
    error_var = "unclassified_diff_se",
    title_y = "Gain in classification (%)",
    breaks = c(0, 10, 20, 30),
    break_labels = c("0", "10", "20", "30"))


candidate_gain_plot <- draw_barplot_simple(
    data = richness_summary,
    x_var = "ecosystem",
    y_var = "mean_richness",
    color = rep("gray", 9),
    error_var = "se_richness",
    title_y = "Candidates richness gain",
    breaks = c(0, 20, 40, 60),
    break_labels = c("0", "20", "40", "60"),
    legend_title = "Ecosystem") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



panel <- plot_grid(
    classification_gain_plot,
    rich_gain_plot,
    candidate_gain_plot,
    align = "v",
    labels = c("A", "B", "C"),
    rel_heights = c(1, 1, 1.3),
    ncol = 1)



#Saave plot
ggsave(
    "results//validation_panel/top_panel.png",
    width = 19,
    height = 18,
    units = "cm")