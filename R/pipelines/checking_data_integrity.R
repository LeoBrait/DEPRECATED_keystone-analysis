# Checking DPAAN and CPR absent in the annotation dataset
library(tidyverse)



if (!dir.exists("results/tables")) {
    dir.create("results/tables")
}


newrad <- radiation %>%
    filter(microgroup != "Bonafide")

all_taxa <- colnames(phyla_abundances[8:length(phyla_abundances)])
all_taxa <- gsub("\\.", " ", all_taxa)

absent_newrad <- newrad %>%
    filter(!taxon %in% all_taxa)

if (nrow(absent_newrad) > 0) {
    print("The following rad-taxa differ in the reference dataset:")
    print(absent_newrad)
}

# Checking samples missing in tha annotation
full_metadata <- read_csv(paste0("data/metadata/", metadata_table))

missing_data <- full_metadata %>%
    filter(!samples %in% phyla_abundances$samples)

write.csv(
    missing_data,
    paste0(
        "results/", analysis_frame, "/missing_data.csv"),
    row.names = FALSE)

write.csv(
    absent_newrad,
    paste0(
        "results/", analysis_frame, "/absent_newrad.csv"),
    row.names = FALSE)

print(paste0(
    "The samples differing in the annotation and metadata",
    "are in results/<analysis_frame>/tables/missing_data.csv"))