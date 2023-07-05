#' @title Calculate cosine similarity between various matrices
#' @description 2 nested functions.
#' First function calculates cosine similarity between 2 matrices.
#' Second function calculates cosine similarity between all matrices in a list.
#' @param matrix1 First matrix
#' @param matrix2 Second matrix
#' @param matrices List of matrices
#' @return Similarity matrix
#' @usage In mostly cases, you will only need to use the second function
#' with a given list of matrices.
#' @require matrixStats
#' @author Bright Mage
library(matrixStats)

cosine_similarity <- function(matrix1, matrix2) {
  flattened1 <- c(matrix1)
  flattened2 <- c(matrix2)
  dot_product <- sum(flattened1 * flattened2)
  norm_product <- sqrt(sum(flattened1^2)) * sqrt(sum(flattened2^2))
  similarity <- dot_product / norm_product
  return(similarity)
}

compute_matrices_similarity <- function(matrices) {
  num_matrices <- length(matrices)
  similarity_matrix <- matrix(0, nrow = num_matrices, ncol = num_matrices)

  for (i in 1:num_matrices) {
    for (j in i:num_matrices) {
      similarity <- cosine_similarity(matrices[[i]], matrices[[j]])
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }

  return(similarity_matrix)
}