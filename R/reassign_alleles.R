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

  # call the TIgGER helper (keeps same arg names as your original)
  db <- reassignAllelesVDJ(db,
    germline_db = seg_references,
    call_column = call_col,
    seq = seq_col,
    trim_seq = trim_seq,
    treat_multigene_as_uncalled = treat_multigene_as_uncalled,
    ignored_regex = ignored_regex
  )

  # If reassignAlleles created the genotype column, count changes and copy back
  if (paste0(seg, "_call_reassigned") %in% names(db)) {
    num_changes <- sum(db[[call_col]] != db[[paste0(seg, "_call_reassigned")]], na.rm = TRUE)
    if (num_changes > 0 && !quiet) {
      message(sprintf("Reassigned %d sequences for %s segment.", num_changes, toupper(seg)))
    }
    db[[call_col]] <- db[[paste0(seg, "_call_reassigned")]]
    db[[paste0(seg, "_call_reassigned")]] <- NULL
  } else {
    warning(sprintf("Expected column '%s_call_reassigned' not found after reassignAlleles for segment '%s'.", seg, toupper(seg)))
  }

  return(db)
}

# Internal helper: compute Hamming distance matrix between query sequences and a list
# of reference sequences. Returns a numeric matrix with rows = queries and
# columns = references.
.compute_hamming_dist_mat <- function(sequences, ref_list, trim_seq = FALSE, data = NULL, germline_trim_columns = NULL, indices = NULL, ignored_regex = "[\\.N-]") {
  # sequences: character vector of full-length sequences (subset already taken)
  # ref_list: named list/vector of reference sequences
  n_q <- length(sequences)
  n_r <- length(ref_list)
  if (n_q == 0 || n_r == 0) {
    return(matrix(numeric(0), nrow = n_q, ncol = n_r))
  }

  dists <- lapply(ref_list, function(x) {
    ref_seqs <- if (trim_seq && !is.null(indices) && !is.null(germline_trim_columns) && !is.null(data)) {
      # repeat and trim the reference sequence for each query according to data
      substr(rep(x, n_q), data[[germline_trim_columns[1]]][indices], data[[germline_trim_columns[2]]][indices])
    } else {
      rep(x, n_q)
    }
    sapply(
      getMutatedPositions(
        sequences,
        ref_seqs,
        match_instead = FALSE,
        ignored_regex = ignored_regex
      ),
      length
    )
  })
  matrix(unlist(dists), ncol = n_r, byrow = FALSE)
}


# Internal helper: for each row of a distance matrix, return indices of minima
.best_match_indices <- function(dist_mat) {
  if (length(dist_mat) == 0) {
    return(list())
  }
  lapply(seq_len(nrow(dist_mat)), function(i) which(dist_mat[i, ] == min(dist_mat[i, ])))
}

#' Reassign alleles by simple Hamming distance to a germline database
#'
#' This function attempts to reassign allele calls (V/D/J) based on the
#' provided `germline_db` by computing simple Hamming distances between the
#' (optionally trimmed) query sequences and candidate germline sequences.
#'
#' @param data A data.frame/tibble containing sequence and call columns.
#' @param germline_db A named list of germline alleles (names are allele ids).
#' @param call_column Character. Column holding original allele calls (e.g. "v_call").
#' @param seq Character. Column holding the alignment used for distance (default "sequence_alignment").
#' @param method Character. Only "hamming" is currently supported.
#' @param trim_seq Logical. If TRUE, reference sequences are trimmed using
#'   *_germline_start/_germline_end positions from `data`.
#' @param keep_gene Character. One of "gene", "family", or "repertoire".
#' @param ignored_regex Regex for characters to ignore when comparing.
#' @param treat_multigene_as_uncalled Logical. Treat multi-gene calls as uncalled.
#' @param top_k Integer. Number of top alleles to keep. If NULL, keep all alleles. Default is 3.
#' @param top_by Character. One of "alphabetical" or "mutation_count". Default is "alphabetical".
#' @return `data` with a new column like "v_call_reassigned" containing the reassigned calls.
#' @export
reassignAllelesVDJ <- function(
  data,
  germline_db,
  call_column = "v_call",
  seq = "sequence_alignment",
  method = "hamming",
  trim_seq = FALSE,
  keep_gene = c("gene", "family", "repertoire"),
  ignored_regex = "[\\.N-]",
  treat_multigene_as_uncalled = FALSE,
  top_k = 3,
  top_by = "alphabetical"
) {
  keep_gene <- match.arg(keep_gene)
  sequences <- as.character(data[[seq]])

  seg <- substr(call_column, 1, 1)

  # Optional trimming of input sequences based on alignment positions
  # useful for D and J reassignment
  if (trim_seq) {
    germline_trim_columns <- c(
      paste0(seg, "_germline_start"),
      paste0(seg, "_germline_end")
    )
    if (seq == "sequence_alignment") {
      seq_trim_columns <- germline_trim_columns
    } else if (seq == "sequence") {
      seq_trim_columns <- c(
        paste0(seg, "_sequence_start"),
        paste0(seg, "_sequence_end")
      )
    } else {
      stop("seq column for trimming must be either sequence_alignment or sequence")
    }

    missing_cols <- setdiff(seq_trim_columns, colnames(data))
    if (length(missing_cols) > 0) {
      stop(
        "Cannot trim sequences: missing columns ",
        paste(missing_cols, collapse = ", ")
      )
    }

    sequences <- substr(
      sequences,
      data[[seq_trim_columns[1]]],
      data[[seq_trim_columns[2]]]
    )
  }

  # Original calls and container for reassigned calls
  calls <- getAllele(data[[call_column]], first = FALSE, strip_d = FALSE)
  calls_reassign <- rep("", length(calls))

  # Identify valid calls (not empty/NA)
  has_call <- !is.na(calls) & nzchar(calls)

  # Detect multi-gene calls when relevant
  if (keep_gene %in% c("gene", "repertoire") && treat_multigene_as_uncalled) {
    genes_full <- getGene(calls, first = FALSE, strip_d = FALSE)
    genes_split <- strsplit(genes_full, ",")
    is_multigene <- vapply(
      genes_split,
      function(x) length(unique(trimws(x))) > 1,
      logical(1)
    )
  } else {
    is_multigene <- rep(FALSE, length(calls))
  }

  # Grouping variable g (by gene/family/repertoire) and germline labels
  if (keep_gene == "gene") {
    g <- getGene(calls, first = TRUE, strip_d = FALSE)
    germline <- getGene(names(germline_db), strip_d = TRUE)
    names(germline) <- names(germline_db)
  } else if (keep_gene == "family") {
    g <- getFamily(calls, first = TRUE, strip_d = FALSE)
    germline <- getFamily(names(germline_db), strip_d = TRUE)
    names(germline) <- names(germline_db)
  } else if (keep_gene == "repertoire") {
    # same behaviour as original code
    g <- rep(calls, length(calls))
    germline <- rep(calls, length(germline_db))
    names(germline) <- names(germline_db)
  } else {
    stop("Unknown keep_gene value: ", keep_gene)
  }

  # Define heterozygous vs homozygous groups in the genotype
  hetero <- unique(germline[which(duplicated(germline))])
  homo <- germline[!(germline %in% hetero)]
  homo_alleles <- names(homo)
  names(homo_alleles) <- homo

  # Homozygous: direct mapping, but skip multi-gene rows
  homo_calls_i <- which(g %in% homo & !is_multigene & has_call)
  if (length(homo_calls_i) > 0) {
    calls_reassign[homo_calls_i] <- homo_alleles[g[homo_calls_i]]
  }

  # Heterozygous: choose best allele(s) within each gene/group, skip multi-gene rows
  if (length(hetero) > 0) {
    for (het in hetero) {
      ind <- which(g %in% het & !is_multigene & has_call)
      if (length(ind) > 0) {
        het_alleles <- names(germline[which(germline == het)])
        het_seqs <- germline_db[het_alleles]

        if (method == "hamming") {
          dist_mat <- .compute_hamming_dist_mat(
            sequences = sequences[ind],
            ref_list = het_seqs,
            trim_seq = trim_seq,
            data = data,
            germline_trim_columns = germline_trim_columns,
            indices = ind,
            ignored_regex = ignored_regex
          )
        } else {
          stop("Only Hamming distance is currently supported as a method.")
        }

        best_match <- .best_match_indices(dist_mat)
        best_alleles <- lapply(best_match, function(x) het_alleles[x])

        if (!is.null(top_k) && !is.na(top_k) && top_by == "alphabetical") {
          best_alleles <- lapply(best_alleles, function(x) {
            if (length(x) > top_k) head(sort(x), top_k) else x
          })
        }

        calls_reassign[ind] <- unlist(lapply(best_alleles, paste, collapse = ","))
      }
    }
  }

  # Rows already handled via homo/hetero (excluding multi-gene rows)
  hetero_calls_i <- which(g %in% hetero & !is_multigene & has_call)

  # Everything else (including multi-gene rows) is treated as "not called"
  # We only reassign if there was an original call
  not_called <- setdiff(which(has_call), c(homo_calls_i, hetero_calls_i))

  if (length(not_called) > 0) {
    if (method == "hamming") {
      dist_mat <- .compute_hamming_dist_mat(
        sequences = sequences[not_called],
        ref_list = germline_db,
        trim_seq = trim_seq,
        data = data,
        germline_trim_columns = germline_trim_columns,
        indices = not_called,
        ignored_regex = ignored_regex
      )
    } else {
      stop("Only Hamming distance is currently supported as a method.")
    }

    best_match <- .best_match_indices(dist_mat)
    best_alleles <- lapply(best_match, function(x) names(germline_db[x]))

    if (!is.null(top_k) && !is.na(top_k) && top_by == "alphabetical") {
      best_alleles <- lapply(best_alleles, function(x) {
        if (length(x) > top_k) head(sort(x), top_k) else x
      })
    }

    calls_reassign[not_called] <- unlist(lapply(best_alleles, paste, collapse = ","))
  }

  # Informative warning if nothing changed
  if (all(calls_reassign == data[[call_column]])) {
    msg <- "No allele assignment corrections made."
    if (all(g %in% homo) && length(hetero) > 0) {
      keep_opt <- eval(formals(reassignAlleles)$keep_gene)
      i <- match(keep_gene, keep_opt)
      rec_opt <- paste(keep_opt[(i + 1):length(keep_opt)], collapse = ", ")
      msg <- paste(msg, "Consider setting keep_gene to one of:", rec_opt)
    }
    warning(msg)
  }

  # Write result column (e.g. v_call_reassigned / d_call_reassigned / j_call_reassigned)
  data[[paste0(seg, "_call_reassigned")]] <- calls_reassign
  return(data)
}
