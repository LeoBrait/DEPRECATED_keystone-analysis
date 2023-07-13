####################### General Community Structures ###########################

microgroups_abundance <- phyla_abundances_long %>%
  mutate(
    microGroup = factor(case_when(
      Taxon %in% dpann_groups ~ "DPANN",
      Taxon %in% cpr_groups ~ "CPR",
      TRUE ~ "Bonafide"
    ))
  ) %>%
  group_by(
    samples,
    life_style,
    biosphere,
    ecosystem,
    habitat,
    latitude,
    longitude,
    microGroup) %>%
  summarise(abundance = sum(abundance)) %>%
  complete(microGroup) %>%
  replace_na(list(abundance = 0)) %>%
  ungroup()

microgroups_richness <- phyla_abundances_long %>%
  mutate(
    microGroup = factor(case_when(
      Taxon %in% dpann_groups ~ "DPANN",
      Taxon %in% cpr_groups ~ "CPR",
      TRUE ~ "Bonafide"
    ))
  ) %>%
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
  ungroup()

microgroups_prevalence <- microgroups_abundance %>%
  left_join(microgroups_richness, by = c(
    "samples",
    "life_style",
    "biosphere",
    "ecosystem",
    "habitat",
    "latitude",
    "longitude",
    "microGroup"
  ))


microgroups_prevalence_persite <- microgroups_prevalence %>%
  group_by(
    biosphere,
    life_style,
    ecosystem,
    habitat,
    latitude,
    longitude,
    microGroup
  ) %>%
  summarise(
    meanAbu = mean(abundance, na.rm = TRUE),
    meanRich = mean(richness, na.rm = TRUE)
  ) %>%
  ungroup()

####################### life-style General Structures ##########################

lifestyle_microgroups_prevalence <- microgroups_prevalence_persite %>%
  group_by(life_style, microGroup) %>%
  summarise(
    n = n(),
    mean_abu_life = mean(meanAbu),
    sd_abu_life = sd(meanAbu),
    se_abu_life = sd_abu_life / sqrt(n),
    mean_rich_life = mean(meanRich),
    sd_rich_life = sd(meanRich),
    se_rich_life = sd_rich_life / sqrt(n)) %>%
    select(-n) %>%
    ungroup()

####################### ecosystem General Structures ##########################

ecosystem_microgroups_prevalence <- microgroups_prevalence_persite %>%
  group_by(ecosystem, microGroup) %>%
  summarise(
    n = n(),
    mean_abu_life = mean(meanAbu),
    sd_abu_life = sd(meanAbu),
    se_abu_life = sd_abu_life / sqrt(n),
    mean_rich_life = mean(meanRich),
    sd_rich_life = sd(meanRich),
    se_rich_life = sd_rich_life / sqrt(n)) %>%
    select(-n) %>%
    ungroup()

################################################################################
####################### Keystone taxa General ##################################

keystones_mean_richness_microgroup <- keystones_abundance_habitat %>%
  filter(LIASP_isKeystone == 1) %>%
  filter(mean_abundance != 0) %>%
  group_by(life_style, ecosystem, habitat, microGroup) %>%
  summarise(
    mean_rich_keystone = n(),
    sd_rich_keystone = sd(mean_rich_keystone),
    se_rich_keystone = sd_rich_keystone / sqrt(n())
  ) %>%
  complete(microGroup) %>%
  replace_na(list(
    mean_rich_keystone = 0,
    sd_rich_keystone = 0,
    se_rich_keystone = 0)) %>%
  ungroup()

################################################################################
####################### Keystone taxa life-style ###############################


lifestyle_keystones_abundances <- keystones_abundance_habitat %>%
  filter(LIASP_isKeystone == 1) %>%
  group_by(life_style, microGroup) %>%
  dplyr::summarise(
    n = n(),
    mean_abu_keystone = mean(mean_abundance),
    sd_abu_keystone = sd(mean_abundance),
    se_abu_keystone = sd_abu_keystone / sqrt(n())) %>%
  ungroup() %>%
  select(-n)

lifestyle_keystones_richness <- keystones_mean_richness_microgroup %>%
  group_by(life_style, microGroup) %>%
  summarise(
    richness = mean(mean_rich_keystone),
    sd_rich_keystone = sd(mean_rich_keystone),
    se_rich_keystone = sd_rich_keystone / sqrt(n())) %>%
  ungroup()

# merge abundance and richness
lifestyle_microgroups_keystones <- lifestyle_keystones_abundances %>%
  left_join(
    lifestyle_keystones_richness,
    by = c("life_style", "microGroup")) %>%
  ungroup()

################################################################################
####################### Keystone taxa ecosystem ################################

ecosystem_keystones_abundances <- keystones_abundance_habitat %>%
  group_by(ecosystem, microGroup) %>%
  filter(LIASP_isKeystone == 1) %>%
  summarise(
    mean_abu_keystone = mean(mean_abundance),
    sd_abu_keystone = sd(mean_abundance),
    se_abu_keystone = sd_abu_keystone / sqrt(n())) %>%
  ungroup() %>%
  complete(ecosystem, microGroup) %>%
  replace_na(list(
    mean_abu_keystone = 0,
    sd_abu_keystone = 0,
    se_abu_keystone = 0))

ecosystem_keystones_richness <- keystones_mean_richness_microgroup %>%
  group_by(ecosystem, microGroup) %>%
  summarise(
    richness = mean(mean_rich_keystone),
    sd_rich_keystone = sd(mean_rich_keystone),
    se_rich_keystone = sd_rich_keystone / sqrt(n())) %>%
  ungroup() %>%
  complete(ecosystem, microGroup) %>%
  replace_na(list(
    richness = 0,
    sd_rich_keystone = 0,
    se_rich_keystone = 0))

# merge abundance and richness
ecosystem_microgroups_keystones <- ecosystem_keystones_abundances %>%
  left_join(
    ecosystem_keystones_richness,
    by = c("ecosystem", "microGroup")) %>%
  complete(microGroup) %>%
  replace_na(list(
    richness = 0,
    sd_rich_keystone = 0,
    se_rich_keystone = 0)) %>%
  ungroup()