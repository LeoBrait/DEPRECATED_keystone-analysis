### Simper Ambient
coresforanalysis <- 6
permutations <- 4999


ecossystems <- keystones_phyla0$level_2
lifestyle <- keystones_phyla0$level_1
input_table <- keystones_phylanumeric

#Simper Ecossystems
pathof_individual_tables <- "outputs//tables//simper_complete//ecossystem//"
nameof_summirizedtable <- "simper_ecossystems_summary.csv"
nameof_finaltable <- "outputs//simper//simper_ecossystems_summary(adjusted).csv"
category <- ecossystems
workpoint <- "ecossystems"
source(
        "source//analysis//simper//simper.r"
      )

#Simper lifestyle
pathof_individual_tables <- "outputs//tables//simper_complete//lifestyle//"
nameof_summirizedtable <- "simper_lifestyle_summary.csv"
nameof_finaltable <- "outputs//tables//simper_lifestyle_summary(adjusted).csv"
category <- lifestyle
workpoint <- "lifestyle"
source(
        "source//simper.r"
      )

################################################################################
keystones_phyla0 <- read_csv(
                    "inputs//phyla.general.relative.matrix_2019-05-09.CSV"
                  )

keystones_phylanumeric <- keystones_phyla0
keystones_phylanumeric[, 1:4] <- list(NULL)

################################################################################

raw_simper <-  simper(
                  input_table,
                  group = category,
                  parallel = coresforanalysis,
                  permutations = permutations
                      )
save.image(paste("workpoints//", workpoint, "simper_step1.RData", sep = ""))


######SAVING ALL TABLES OF SPECIFIC SIMPER
summary_simper <- summary(raw_simper)
for (name_of_table in names(summary_simper)){
            data <- summary_simper[[name_of_table]]
            write.csv(file = paste(
                                    pathof_individual_tables,
                                    name_of_table,
                                    ".csv",
                                    sep = ""
                                  ),
                       x = data
                      )
                                }

load("work_environments//simper.RData")

#Filtering the data
simpertables <- list()

##testing
df <- summary_simper[["saline water_soil"]]
View(df)



for(name_of_table in names(summary_simper)){
                #ambient setting
                df <- summary_simper[[name_of_table]]
                prev <- 0
                #getting the cumsums and converting to individual simper values
                for(rows in 1:nrow(df)){
                                df[rows, 8] <-  df[rows, 6] - prev
                                prev <-  df[rows, 6]
                                        }

                #setting names for the data.frame and setting comparisons
                df$comparison <- rep(name_of_table, nrow(df))
                colnames(df) <- c(
                        "average", "sd", "ratio", "ava",
                        "avb", "cumsum", "p", "contribution", "comparison"
                                 )
                df <- tibble::rownames_to_column(df, "OTU")
                simpertables[[name_of_table]] <- df
                }



#Reuniting all data in a single data.frame
simper_clean <-  bind_rows(simpertables)
simper_clean <- simper_clean[simper_clean$p < 0.05, ]
View(simper_clean)



write.csv(file = paste(
                        "outputs//simper//",
                        nameof_summirizedtable,
                        sep = ""
                      ),
          x = simper_clean
         )

df0 <- simper_clean
#AMBIENT
df0 <- as.data.frame(cbind(df0$comparison, df0$contribution, df0$OTU))
colnames(df0) <- c("comparison", "contribution", "OTU")
df0$contribution <- as.numeric(df0$contribution)
levels <- unique(df0$comparison)
subset_contribution <- list(NULL)
firsts <- c()
seconds <- c()
thirds <- c()
fourths <- c()
fives <- c()
six <- c()
seven <- c()
eight <- c()
nine <- c()

#listing subset by ist factors and reordering and loading
for (x in 1:length(levels)){
        dummy <- subset(df0, comparison == levels[x])
        dummy <- arrange(dummy, -contribution)
        row.names(dummy) <- NULL
        subset_contribution[x] <- list(dummy)
                          }




for (x in 1:length(subset_contribution)) {
        firsts[length(firsts) + 1] <- as.data.frame(subset_contribution[x])[1, 3]
        print(as.data.frame(subset_contribution[x])[1, 3])
                                                                                      }



for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[2, 3]
        if (dummy %in% firsts) {
                } else {
                        seconds[length(seconds) + 1] <- as.data.frame(subset_contribution[x])[2, 3]
                        print(as.data.frame(subset_contribution[x])[2, 3])
                                                                                                        }
                                                                                                          }

for (x in 1:length(subset_contribution)) {
        dummy <- as.data.frame(subset_contribution[x])[3, 3]
        if (dummy %in% firsts || dummy %in% seconds) {
                } else {
                        thirds[length(thirds) + 1] <- as.data.frame(subset_contribution[x])[3, 3]
                        print(as.data.frame(subset_contribution[x])[3, 3])
                                                                                                      }
                                                                                                        }

for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[4, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds) {
                } else {
                        fourths[length(fourths) + 1] <- as.data.frame(subset_contribution[x])[4, 3]
                        print(as.data.frame(subset_contribution[x])[4, 3])
                                                                                                      }
                                                                                                        }


for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[5, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds || dummy %in% fourths) {
                } else {
                        fives[length(fives) + 1] <- as.data.frame(subset_contribution[x])[5, 3]
                        print(as.data.frame(subset_contribution[x])[5, 3])
                                                                                                      }
                                                                                                        }
for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[6, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds || dummy %in% fourths || dummy %in% fives) {
                } else {
                        six[length(six) + 1] <- as.data.frame(subset_contribution[x])[6, 3]
                        print(as.data.frame(subset_contribution[x])[6, 3])
                                                                                                      }
                                                                                                        }

for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[7, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds || dummy %in% fourths || dummy %in% fives || dummy %in% six) {
                } else {
                        seven[length(seven) + 1] <- as.data.frame(subset_contribution[x])[7, 3]
                        print(as.data.frame(subset_contribution[x])[7, 3])
                                                                                                      }
                                                                                                        }
for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[8, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds || dummy %in% fourths || dummy %in% fives || dummy %in% six || dummy %in% seven) {
                } else {
                        eight[length(eight) + 1] <- as.data.frame(subset_contribution[x])[8, 3]
                        print(as.data.frame(subset_contribution[x])[8, 3])
                                                                                                      }
                                                                                                        }
for (x in 1:length(subset_contribution)){
        dummy <- as.data.frame(subset_contribution[x])[9, 3]
        if (dummy %in% firsts || dummy %in% seconds || dummy %in% thirds || dummy %in% fourths || dummy %in% fives || dummy %in% six || dummy %in% seven || dummy %in% eight) {
                } else {
                        nine[length( nine) + 1] <- as.data.frame(subset_contribution[x])[9, 3]
                        print(as.data.frame(subset_contribution[x])[9, 3])
                                                                                                      }
                                                                                                        }


firsts <- as.data.frame(table(firsts))
colnames(firsts) <- c("OTU", "Freq")
firsts <- arrange(firsts, -Freq)
firsts$pos <- rep(1, nrow(firsts))

seconds <- as.data.frame(table(seconds))
colnames(seconds) <- c("OTU", "Freq")
seconds <- arrange(seconds, -Freq)
seconds$pos <- rep(2, nrow(seconds))

thirds <- as.data.frame(table(thirds))
colnames(thirds) <- c("OTU", "Freq")
thirds <- arrange(thirds, -Freq)
thirds$pos <- rep(3, nrow(thirds))

fourths <- as.data.frame(table(fourths))
colnames(fourths) <- c("OTU", "Freq")
fourths <- arrange(fourths, -Freq)
fourths$pos <- rep(4, nrow(fourths))

fives <- as.data.frame(table(fives))
colnames(fives) <- c("OTU", "Freq")
fives <- arrange(fives, -Freq)
fives$pos <- rep(5, nrow(fives))

six <- as.data.frame(table(six))
colnames(six) <- c("OTU", "Freq")
six  <- arrange(six, -Freq)
six$pos <- rep(6, nrow(six))

seven <- as.data.frame(table(seven))
colnames(seven) <- c("OTU", "Freq")
seven  <- arrange(seven, -Freq)
seven$pos <- rep(7, nrow(seven))

eight <- as.data.frame(table(eight))
colnames(eight) <- c("OTU", "Freq")
eight  <- arrange(eight, -Freq)
eight$pos <- rep(8, nrow(eight))

nine <- as.data.frame(table(nine))
colnames(nine) <- c("OTU", "Freq")
nine  <- arrange(nine, -Freq)
nine$pos <- rep(9, nrow(nine))

ranking <- data.table::rbindlist(list(firsts, seconds, thirds, fourths, fives, six, seven, eight, nine))
row.names(ranking) <- NULL
rm(firsts, seconds, thirds, fourths, levels, subset_contribution, df0)
View(ranking)
write.csv(ranking, "results//ranking_simper.csv")