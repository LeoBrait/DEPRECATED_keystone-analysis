#' merge annotation and metadata
#' Select the respective and desired metadata to include into annotation.
#' You will need a "samples" column in both tables.
#'@title merge annotation metadata
#'@param annotation_df the annotated table, with a "samples" column.
#'@param metadata_df the metadata table, with a "samples" column.
#'@param metadata_variables the variables to include in the annotation.
#'@returns The merged table, with the selected metadata variables.
merge_annotation_metadata <- function(
    annotation_df,
    metadata_df,
    metadata_variables) {

    metadata_filtered <- metadata_df %>%
        dplyr::select(
            metadata_variables)

    metagenomes_metadata <- inner_join(
        metadata_filtered,
        annotation_df,
        by = "samples")


return(metagenomes_metadata)
}