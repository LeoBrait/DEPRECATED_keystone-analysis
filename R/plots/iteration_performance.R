library("tidyverse")
library("ggplot2")

data <- read_csv(
    "data/summaries/performance_fastspar_iterations/animal_host-associated.aqueous_humour.csv")
df_long <- gather(data, key = "column_name", value = "values")


ggplot(df_long, aes(x = column_name, y = values)) +
  geom_boxplot() +
  labs(title = "Boxplot of Values", x = "Column Name", y = "Values")



ggplot(df_long, aes(x = column_name, y = values)) +
  geom_point() +
  labs(title = "Scatterplot of Values", x = "Column Name", y = "Values")
