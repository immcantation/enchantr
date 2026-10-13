# The tigger_bayesian_genotype report as it was before the unified genotype
# report, kept for callers that use it (e.g. nf-core/airrflow).

#' Create a TIGGER Bayesian genotype project template
#'
#' Create the directory skeleton for a TIGGER Bayesian genotype project from
#' the package's bundled RStudio project template. This is a convenience
#' wrapper used by RStudio's New Project > New Directory > From Package
#' templates.
#'
#' @param path Path to the directory where the project will be created.
#' @param ... Additional arguments (currently unused).
#' @export
tigger_bayesian_genotype_project <- function(path, ...) {
  skeleton_dir <- file.path(
    system.file(package = "enchantr"), "rstudio",
    "templates", "project",
    "tigger_bayesian_genotype_project_files"
  )
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}

#' Genotype inference used by the tigger_bayesian_genotype report
#'
#' The infer_genotype() of enchantr before the unified genotype report, kept so
#' the tigger_bayesian_genotype report produces the same output.
#' @noRd
.infer_genotype_bayesian <- function(
  db,
  seg,
  loci,
  genotypeby = NULL,
  references = NULL,
  single_assignments = FALSE,
  find_unmutated = FALSE,
  seq = "sequence_alignment",
  priors = c(0.6, 0.4, 0.4, 0.35, 0.25, 0.25, 0.25, 0.25, 0.25)
) {
  call_col <- paste0(seg, "_call")
  seg_references <- NULL
  if (!is.null(references) && find_unmutated == TRUE) {
    # collect reference alleles for this segment across loci
    seg_references <- unlist(lapply(loci, function(locus) {
      # be tolerant if the locus or segment entry is missing
      if (!is.list(references[[locus]]) || is.null(references[[locus]][[toupper(seg)]])) {
        return(NA_character_)
      }
      references[[locus]][[toupper(seg)]]
    }))
    seg_references <- seg_references[!is.na(seg_references) & nzchar(seg_references)]
  }

  if (single_assignments) {
    db <- db[!grepl(",", db[[call_col]]), ]
  }

  # Remove NA values
  db <- db[!is.na(db[[call_col]]), ]

  # Run Bayesian inference on one repertoire (or sub-group); NULL on failure.
  run_one <- function(d) {
    g <- try(
      tigger::inferGenotypeBayesian(
        d,
        find_unmutated = find_unmutated,
        germline_db = seg_references,
        v_call = call_col,
        seq = seq,
        priors = priors,
        genotyped_alleles = TRUE
      )
    )
    if (inherits(g, "try-error")) NULL else g
  }

  # Whole-repertoire inference
  if (is.null(genotypeby)) {
    return(run_one(db))
  }

  # Per-group inference: tag each result with the grouping value and bind
  parts <- lapply(unique(db[[genotypeby]]), function(val) {
    g <- run_one(db[db[[genotypeby]] == val, ])
    if (is.data.frame(g)) {
      g[[genotypeby]] <- val
    }
    g
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) > 0) do.call(rbind, parts) else NULL
}

#' Summarize single-assignment support for genotype genes
#'
#' Internal helper that counts raw productive sequences and unique clones
#' supporting each genotyped gene. Counts are restricted to single-assignment
#' calls so they can be compared against TIgGER's ambiguity-weighted `counts`
#' column.
#'
#' @param db A data.frame/tibble of productive repertoire records.
#' @param genotypeby Optional grouping column name used for per-group genotype
#'   inference.
#' @param clone_id_column Optional clone identifier column name.
#' @return A data.frame with `gene`, `allele`, `sequence_count`,
#'   `clone_count`, and the grouping column when `genotypeby` is provided.
#' @keywords internal
.single_assignment_support_counts <- function(db, genotypeby = NULL, clone_id_column = NULL) {
  if (is.null(db) || nrow(db) == 0) {
    return(data.frame())
  }

  count_by_segment <- lapply(c("v", "d", "j"), function(seg) {
    call_col <- paste0(seg, "_call")
    if (!call_col %in% colnames(db)) {
      return(NULL)
    }

    seg_db <- db[!is.na(db[[call_col]]) & nzchar(db[[call_col]]) & !grepl(",", db[[call_col]]), , drop = FALSE]
    if (nrow(seg_db) == 0) {
      return(NULL)
    }

    seg_db$gene <- alakazam::getGene(seg_db[[call_col]], strip_d = FALSE)
    seg_db$allele <- sub(".*\\*", "", seg_db[[call_col]])
    group_cols <- c(if (!is.null(genotypeby) && genotypeby %in% colnames(seg_db)) genotypeby, "gene", "allele")

    sequence_counts <- seg_db %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::summarise(sequence_count = dplyr::n(), .groups = "drop")

    if (!is.null(clone_id_column) &&
      nzchar(clone_id_column) &&
      clone_id_column %in% colnames(seg_db)) {
      clone_counts <- seg_db %>%
        dplyr::filter(!is.na(.data[[clone_id_column]]), .data[[clone_id_column]] != "") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(clone_count = dplyr::n_distinct(.data[[clone_id_column]]), .groups = "drop")

      dplyr::left_join(sequence_counts, clone_counts, by = group_cols)
    } else {
      sequence_counts$clone_count <- NA_integer_
      sequence_counts
    }
  })

  dplyr::bind_rows(count_by_segment)
}

#' Add interpretation-oriented columns to genotype output
#'
#' Internal helper used by the genotype report to rename display columns,
#' attach raw single-assignment support counts, and enforce a stable
#' presentation order for report tables and TSV outputs.
#'
#' @param genotypes A data.frame returned by `infer_genotype()`.
#' @param db_support A productive repertoire data.frame used to compute support
#'   counts before clone-representative downsampling.
#' @param genotypeby Optional grouping column name.
#' @param clone_id_column Optional clone identifier column name.
#' @return A data.frame with added `sequence_count`/`clone_count`, renamed
#'   `candidate_alleles`, and reordered columns.
#' @keywords internal
.augment_genotype_table <- function(genotypes, db_support, genotypeby = NULL, clone_id_column = NULL) {
  if (is.null(genotypes) || nrow(genotypes) == 0) {
    return(genotypes)
  }

  if ("alleles" %in% colnames(genotypes) && !"candidate_alleles" %in% colnames(genotypes)) {
    genotypes <- dplyr::rename(genotypes, candidate_alleles = alleles)
  }

  support_counts <- .single_assignment_support_counts(
    db = db_support,
    genotypeby = genotypeby,
    clone_id_column = clone_id_column
  )
  count_group_cols <- c(if (!is.null(genotypeby) && genotypeby %in% colnames(genotypes)) genotypeby, "gene")
  has_clone_counts <- nrow(support_counts) > 0 && !all(is.na(support_counts$clone_count))

  summarize_counts <- function(row_idx, count_col) {
    allele_string <- genotypes[row_idx, "genotyped_alleles"]
    if (is.na(allele_string) || allele_string == "") {
      return(NA_character_)
    }

    row_counts <- support_counts
    for (col in count_group_cols) {
      row_counts <- row_counts[row_counts[[col]] == genotypes[row_idx, col], , drop = FALSE]
    }

    alleles <- trimws(unlist(strsplit(allele_string, ",")))
    values <- vapply(alleles, function(allele) {
      hit <- row_counts[row_counts$allele == allele, count_col, drop = TRUE]
      if (length(hit) == 0 || all(is.na(hit))) {
        if (identical(count_col, "clone_count") && !has_clone_counts) {
          return(NA_character_)
        }
        return("0")
      }
      as.character(hit[[1]])
    }, character(1))

    if (identical(count_col, "clone_count") && !has_clone_counts) {
      return(NA_character_)
    }
    paste(values, collapse = ",")
  }

  genotypes$sequence_count <- vapply(seq_len(nrow(genotypes)), summarize_counts, character(1), count_col = "sequence_count")
  genotypes$clone_count <- vapply(seq_len(nrow(genotypes)), summarize_counts, character(1), count_col = "clone_count")

  ordered_cols <- c(
    if (!is.null(genotypeby) && genotypeby %in% colnames(genotypes)) genotypeby,
    "gene",
    "genotyped_alleles",
    "sequence_count",
    "clone_count",
    "k_diff",
    "note",
    "candidate_alleles",
    "counts",
    "total",
    "kh",
    "kd",
    "kt",
    "kq"
  )
  ordered_cols <- intersect(ordered_cols, colnames(genotypes))
  remaining_cols <- setdiff(colnames(genotypes), ordered_cols)

  genotypes[, c(ordered_cols, remaining_cols), drop = FALSE]
}
