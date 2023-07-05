#'@title convert_matrix_tonumeric
#'@description Converts a matrix of chr to numeric
convert_matrix_tonumeric <- function(df) {
  char_cols <- sapply(df, is.character)
  df[char_cols] <- lapply(df[char_cols], as.numeric)
  df
}