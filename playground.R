library("tidyverse")
options(scipen = 9999999)
options(digits = 20)


taxons <- read_csv("data/taxon_abundances/kraken_custom_phyla.csv")
taxons[1] <- NULL

min_nonzero <- function(x) {
  min_val <- min(x[x != 0])
  if (is.finite(min_val)) {
    return(min_val)
  } else {
    return(NA)
  }
}

# Find the minimum non-zero values in each column
min_values_nonzero <- apply(taxons, 2, min_nonzero)

# Print the minimum non-zero values
print(min_values_nonzero)
min_values_nonzero <- as.data.frame(min_values_nonzero)

#transform into csv
write_csv(min_values_nonzero, "min_nonzero.csv")
