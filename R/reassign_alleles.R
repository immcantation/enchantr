#' Create an reassign alleles project
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  reassign alleles
#' to create the skeleton of an reassign alleles project
#' @param  path path to the directory where the project will be created
reassign_alleles_project <- function(path, ...) {
  skeleton_dir <- file.path(
    system.file(package = "enchantr"), "rstudio",
    "templates", "project",
    "reassign_alleles_project_files"
  )
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}

#' Reassign alleles for a single immunoglobulin segment
#'
#' Collects germline alleles for a given segment (V/D/J) across the provided
#' loci, calls TIgGER's `reassignAlleles()` to compute genotyped assignments,
#' copies the genotyped column back to the original `<seg>_call` column and
#' removes the temporary genotype column.
#'
#' @param db A data.frame or tibble with repertoire annotations. Must contain
#'   a column named like "<seg>_call" (for example, "v_call").
#' @param seg Character scalar. Segment name to reassign; typically "v", "d",
#'   or "j".
#' @param loci Character vector of locus names to extract alleles from in
#'   `references`.
#' @param references A nested list of reference alleles indexed by locus. Each
#'   locus should contain elements named "V", "D" or "J" (upper-case).
#' @param treat_multigene_as_uncalled Logical scalar. If `TRUE`, sequences with
#'   multi-gene calls are treated as uncalled during reassignment. Default is
#'   `TRUE`.
#' @param quiet Logical scalar. If `TRUE`, suppress messages about the number of
#'   reassigned sequences. Default is `FALSE`.
#' @return The input `db` with the `<seg>_call` column updated when
#'   genotyped assignments are present; the temporary `<seg>_call_genotyped`
#'   column is removed.
#' @details The function is tolerant of missing loci / missing segment entries
#' in `references` and will warn and skip the segment if no alleles are found.
#' If `reassignAlleles()` does not create the expected genotype column a
#' warning is emitted and `db` is returned unchanged for that segment.
#' @export
reassign_segment <- function(db, seg, loci, references, treat_multigene_as_uncalled = TRUE, quiet = FALSE) {
  # collect reference alleles for this segment across loci
  seg_references <- unlist(lapply(loci, function(locus) {
    # be tolerant if the locus or segment entry is missing
    if (!is.list(references[[locus]]) || is.null(references[[locus]][[toupper(seg)]])) {
      return(NA_character_)
    }
    references[[locus]][[toupper(seg)]]
  }))
  seg_references <- seg_references[!is.na(seg_references) & nzchar(seg_references)]

  if (length(seg_references) == 0) {
    warning(sprintf("No reference alleles found for segment '%s' (loci: %s) — skipping.", seg, paste(loci, collapse = ", ")))
    return(db)
  }

  call_col <- paste0(seg, "_call")
  trim_seq <- seg != "v"
  seq_col <- "sequence"
  ignored_regex <- TRUE
  if (seg == "v" && "sequence_alignment" %in% colnames(db)) {
    seq_col <- "sequence_alignment"
    ignored_regex <- "[\\.N-]"
  }

  original_calls <- db[[call_col]]
  db <- tigger::reassignAlleles(db,
    genotype_db = seg_references,
    v_call = call_col,
    seq = seq_col,
    trim_seq = trim_seq,
    overwrite = TRUE,
    treat_multigene_as_uncalled = treat_multigene_as_uncalled,
    ignored_regex = ignored_regex,
    top_k = 3,
    top_by = "alphabetical",
    strip_d = FALSE,
    reassign_uncalled = FALSE
  )

  num_changes <- sum(original_calls != db[[call_col]], na.rm = TRUE)
  if (num_changes > 0 && !quiet) {
    message(sprintf("Reassigned %d sequences for %s segment.", num_changes, toupper(seg)))
  }

  return(db)
}

#' Reassign alleles per locus and recombine
#'
#' Reassigns each locus' sequences against only that locus' alleles (so a gene
#' absent from the genotype is never compared across loci), running only the
#' segments applicable to the locus' chain (no D for light chains), then
#' row-binds the loci back together in the original row order.
#'
#' @param db AIRR \code{data.frame} with a \code{locus} column.
#' @param loci Character vector of locus names present in \code{db}.
#' @param segments Character vector of segments to reassign (subset of
#'   \code{c("v", "d", "j")}).
#' @param references Nested list of reference alleles indexed by locus and
#'   segment (as returned by \code{dowser::readIMGT}).
#' @return \code{db} with reassigned \code{<seg>_call} columns, in the original row order.
#' @export
reassign_loci <- function(db, loci, segments, references) {
  db[[".ord"]] <- seq_len(nrow(db))
  db <- dplyr::bind_rows(lapply(loci, function(locus) {
    locus_db <- db[db[["locus"]] == locus, ]
    locus_segments <- intersect(
      segments,
      if (isHeavyChain(locus)) c("v", "d", "j") else c("v", "j")
    )
    for (seg in locus_segments) {
      locus_db <- reassign_segment(locus_db, seg, locus, references)
    }
    locus_db
  }))
  db <- db[order(db[[".ord"]]), ]
  db[[".ord"]] <- NULL
  db
}
