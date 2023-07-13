libs_path <- paste0("test_libs/")
if (!file.exists(libs_path)) {
  dir.create(libs_path)
}


source("R/src/install_and_load.R")
install_and_load(
  libs = c(
    "tidyverse" = "2.0.0",
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
    "scales" = "1.2.1"),
    loc = libs_path)