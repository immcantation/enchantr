#' Select clonal representatives with minimal mutations
#'
#' This function selects a representative sequence for each clone, defined as the sequence with the fewest mutations compared to the germline. It supports custom column names for sequence, germline, and clone ID, and adds a flag for the representative.
#'
#' @param data A data.frame containing the repertoire data.
#' @param sequence_col Name of the column with the sequence alignment (default: 'sequence_alignment').
#' @param germline_col Name of the column with the germline alignment (default: 'germline_alignment_d_mask').
#' @param clone_id_col Name of the column with the clone ID (default: 'clone_id').
#' @param parallel Logical, whether to use parallel computation (default: TRUE).
#' @param return_count Logical, whether to return mutation counts (default: TRUE).
#' @return The input data.frame with a new column 'mut' (mutation count) and a logical 'clone_representative' flag.
#' @examples
#' \dontrun{
#' df <- select_clonal_representatives(df)
#' }
#' @export
select_clonal_representatives <- function(
  data,
  sequence_col = "sequence_alignment",
  germline_col = "germline_alignment_d_mask",
  clone_id_col = "clone_id",
  parallel = TRUE,
  return_count = TRUE
) {
  stopifnot(sequence_col %in% names(data))
  stopifnot(germline_col %in% names(data))
  stopifnot(clone_id_col %in% names(data))

  data <- data %>%
    mutate(
      mut = mutation_count(
        .data[[sequence_col]],
        .data[[germline_col]],
        parallel = parallel,
        return_count = return_count
      )
    ) %>%
    group_by(.data[[clone_id_col]]) %>%
    mutate(
      clone_size = n(),
      clone_representative = (.data$mut == min(.data$mut, na.rm = TRUE))
    ) %>%
    ungroup()

  data
}
