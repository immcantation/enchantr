#' Select clonal representatives with minimal mutations
#'
#' This function selects a representative sequence for each clone, defined as the sequence with the fewest mutations compared to the germline. It supports custom column names for sequence, germline, and clone ID, and adds a flag for the representative.
#'
#' @param data A data.frame containing the repertoire data.
#' @param sequence_col Name of the column with the sequence alignment (default: 'sequence_alignment').
#' @param germline_col Name of the column with the germline alignment (default: 'germline_alignment_d_mask').
#' @param clone_id_col Name of the column with the clone ID (default: 'clone_id').
#' @details A sequence with no clone ID is its own representative.
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
  clone_id_col = "clone_id"
) {
  stopifnot(sequence_col %in% names(data))
  stopifnot(germline_col %in% names(data))
  stopifnot(clone_id_col %in% names(data))

  data <- data %>%
    mutate(
      mut = alakazam::seqMismatchCount(
        .data[[sequence_col]],
        .data[[germline_col]]
      ),
      unclonal = is.na(.data[[clone_id_col]]) |
        !nzchar(as.character(.data[[clone_id_col]]))
    ) %>%
    group_by(.data[[clone_id_col]]) %>%
    mutate(
      clone_size = if (any(.data$unclonal)) 1L else n(),
      clone_representative = if (any(.data$unclonal)) TRUE
                             else .data$mut == min(.data$mut, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(-"unclonal")

  data
}
