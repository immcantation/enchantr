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

# Helper function to get most likely genotyped alleles
#' Helper: extract most likely genotyped alleles from inferGenotypeBayesian output
#'
#' Internal helper used to convert columns returned by
#' `tigger::inferGenotypeBayesian()` into a comma-separated list of the
#' most-likely alleles.
#'
#' @param row A single-row character/numeric vector where the first element is
#'   a comma-separated allele list and elements 2:5 are likelihood scores.
#' @return A character string with the top N alleles (comma-separated).
#' @keywords internal
.genotyped_alleles <- function(row) {
  m <- which.max(as.numeric(row[2:5]))
  paste0(unlist(strsplit(row[1], ","))[1:m], collapse = ",")
}

## wrapper to infer genotype by segment, and by genotypeby column and return the results
#' Infer genotype (Bayesian) for a segment
#'
#' Wrapper around `tigger::inferGenotypeBayesian()` that runs genotype
#' inference for a single segment (V/D/J). Optionally, the inference can be
#' run per-group by providing `genotypeby`, in which case a named list of
#' genotype data.frames is returned (one per group).
#'
#' @param db A data.frame or tibble containing rearrangement records (AIRR format).
#' @param seg Character scalar. Segment short name, e.g. "v", "d", or "j".
#' @param loci Character vector. Loci to extract reference alleles from when
#'   `references` is provided (e.g. c("IGH")).
#' @param genotypeby Optional character column name in `db` to run inference per-group.
#' @param references Optional nested list of reference sequences (as produced by dowser/readIMGT). If provided,
#'   alleles for the selected `seg` are collected from each `locus` in `loci` and
#'   passed to `tigger::inferGenotypeBayesian()` via `germline_db`.
#' @param find_unmutated Logical. Passed to `inferGenotypeBayesian()`.
#' @param single_assignments Logical. Passed to `inferGenotypeBayesian()`. If TRUE, only sequences with single allele assignments (no commas) are used for inference.
#' @param seq Column name holding alignment sequences (default: "sequence_alignment").
#' @param priors Numeric vector of prior probabilities for genotype model; see
#'   `tigger::inferGenotypeBayesian()` for details.
#' @return If `genotypeby` is NULL, returns a single data.frame (the genotype table).
#' If `genotypeby` is provided, returns a named list of genotype tables (one per group).
#' @export
infer_genotype <- function(
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

  if (!is.null(genotypeby)) {
    # run inference per-group, collect genotype tables and per-group references
    groups <- unique(db[[genotypeby]])
    genotype_list <- list()
    for (val in groups) {
      db_genotypeby <- db[db[[genotypeby]] == val, ]
      genotype_try <- try(
        tigger::inferGenotypeBayesian(
          db_genotypeby,
          find_unmutated = find_unmutated,
          germline_db = seg_references,
          v_call = call_col,
          seq = seq,
          priors = priors
        )
      )
      if (inherits(genotype_try, "try-error")) {
        genotype_list[[as.character(val)]] <- NULL
        next
      }
      genotype_try$genotyped_alleles <- apply(genotype_try[, c(2, 6:9)], 1, .genotyped_alleles)
      gt_df <- genotype_try
      if (is.data.frame(gt_df)) {
        gt_df[[genotypeby]] <- val
      }
      genotype_list[[as.character(val)]] <- gt_df
    }
    genotype_bound <- NULL
    non_null_gts <- genotype_list[!vapply(genotype_list, is.null, logical(1))]
    if (length(non_null_gts) > 0) {
      genotype_bound <- do.call(rbind, non_null_gts)
    }
    genotype <- genotype_bound
  } else {
    genotype <- try(
      tigger::inferGenotypeBayesian(
        db,
        find_unmutated = find_unmutated,
        germline_db = seg_references,
        v_call = call_col,
        seq = seq,
        priors = priors
      )
    )

    if (inherits(genotype, "try-error")) {
      return(NULL)
    }

    genotype$genotyped_alleles <- apply(genotype[, c(2, 6:9)], 1, .genotyped_alleles)
  }

  return(genotype)
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

#' Write genotype report(s) to disk
#'
#' Write genotype result tables to `out_dir`. If `genotypeby` is provided,
#' write one file per group named `<group>_<out_fn>`, otherwise write `out_fn`.
#'
#' @param genotypes Data.frame of genotype results.
#' @param out_dir Directory to write files into.
#' @param out_fn Filename to use when `genotypeby` is NULL, otherwise used as suffix.
#' @param genotypeby Optional column name in `genotypes` that contains grouping values.
#' @export
write_genotypes <- function(genotypes, out_dir, out_fn, genotypeby = NULL) {
  if (!is.null(genotypeby)) {
    groups <- unique(genotypes[[genotypeby]])
    for (val in groups) {
      genotypes_genotypeby <- genotypes[genotypes[[genotypeby]] == val, ]
      write.table(
        genotypes_genotypeby,
        file = file.path(out_dir, paste0(val, "_", out_fn)),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
    }
  } else {
    write.table(
      genotypes,
      file = file.path(out_dir, out_fn),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }
}

# Helper to build personal germline reference from genotype table
#' Helper: build personal germline reference from genotype table
#'
#' Internal helper used to construct a personal germline reference sequence
#' list from a genotype table and the original reference sequences.
#'' @param genotype A data.frame with genotype results for a segment.
#' @param reference Named character/list of original germline sequences (names like "GENE*ALLELE").
#' @return A named character/list of germline sequences for the personal genotype.
#' @keywords internal
.genotyped_reference <- function(genotype, reference) {
  if (is.null(reference)) {
    return(NULL)
  }
  not_in_genotype <- !(sapply(strsplit(names(reference), "*", fixed = TRUE), `[`, 1) %in% genotype$gene)
  genotype_reference <- reference[not_in_genotype]
  # Add genotyped alleles
  for (i in seq_len(nrow(genotype))) {
    gene <- genotype[i, "gene"]
    candidate_col <- if ("candidate_alleles" %in% colnames(genotype)) "candidate_alleles" else "alleles"
    alleles <- if (genotype[i, "genotyped_alleles"] == "") genotype[i, candidate_col] else genotype[i, "genotyped_alleles"]
    alleles <- unlist(strsplit(alleles, ","))
    ind <- names(reference) %in% paste(gene, alleles, sep = "*")
    genotype_reference <- c(genotype_reference, reference[ind])
  }

  genotype_reference
}

#' Generate genotype-informed reference FASTA files
#'
#' For a combined genotype table (optionally containing a grouping column
#' given by `genotypeby`), build per-locus, per-segment genotype-aware
#' reference FASTA files. When `genotypeby` is provided, a separate
#' output tree is created for each group under `file.path(output_dir, group, species, "vdj")`.
#'
#' The function copies the baseline files from `references_dir` into the
#' target location and overwrites locus FASTA files with genotype-filtered
#' sequences derived from `genotypes` and `references`.
#'
#' @param genotypes A data.frame with genotype results. Must contain a column `gene` and,
#'   when `genotypeby` is set, a grouping column with that name.
#' @param references Nested list of original reference sequences (as returned by
#'   dowser::readIMGT), indexed by locus and segment (e.g. `references[["IGH"]][["V"]]`).
#' @param output_dir Root output directory where per-group reference directories will be created.
#' @param references_dir Path to the original IMGT reference directory to copy baseline files from (the function
#'   will copy files found under this directory into each group's `.../vdj` folder).
#' @param loci Character vector of locus names to process (e.g. `c("IGH")`).
#' @param species Species name used to place files under `<group>/<species>/vdj` (default: "human").
#' @param genotypeby Optional string naming the grouping column in `genotypes`. If NULL a single reference is written.
#' @param quiet Logical; if TRUE minimize console messages.
#' @return Invisibly returns NULL. The primary effect is to write FASTA files to disk.
#' @export
generate_genotyped_reference <- function(
  genotypes,
  references,
  output_dir,
  germline_dir,
  references_dir,
  loci,
  species,
  genotypeby = NULL,
  quiet = FALSE
) {
  if (!is.null(genotypeby)) {
    groups <- unique(genotypes[[genotypeby]])
    for (val in groups) {
      genotypes_genotypeby <- genotypes[genotypes[[genotypeby]] == val, ]

      # create output dir for this group
      output_dir_genotypeby <- file.path(output_dir, as.character(val), germline_dir, species, "vdj")
      dir.create(output_dir_genotypeby, showWarnings = FALSE, recursive = TRUE)

      # copy the files from the original references dir
      file.copy(list.files(references_dir, full.names = TRUE), output_dir_genotypeby, recursive = TRUE)

      # get the loci files
      loci_files <- list.files(output_dir_genotypeby, paste0(loci, collapse = "|"), full.names = TRUE)

      # for each locus, and segment, build the reference and save it
      for (locus in loci) {
        locus_file <- loci_files[grep(locus, loci_files)]
        if (length(locus_file) == 0) {
          next
        }
        # get the locus genotype
        genotypes_genotypeby_locus <- genotypes_genotypeby[grepl(locus, genotypes_genotypeby$gene), ]

        for (seg in c("V", "D", "J")) {
          call <- paste0(locus, seg)
          # check that there are genes for the locus and segment
          if (any(grepl(call, genotypes_genotypeby_locus$gene))) {
            seg_references <- references[[locus]][[seg]]
            # get the genotype reference
            seg_ref_group <- .genotyped_reference(genotypes_genotypeby_locus[grepl(call, genotypes_genotypeby_locus$gene), ], seg_references)
            # write the genotype reference
            seg_ref_file <- locus_file[grepl(paste0(locus, seg, "\\."), basename(locus_file))]
            if (length(seg_ref_file) != 1) {
              stop(sprintf("Expected exactly one reference file for %s, found %d.", call, length(seg_ref_file)))
            }
            tigger::writeFasta(seg_ref_group, file.path(output_dir_genotypeby, basename(seg_ref_file)))
            if (!quiet) {
              message(sprintf("%s Genotyped Reference database for %s: [%s](%s)\n", call, val, seg_ref_file, seg_ref_file))
            }
          } else {
            next
          }
        }
      }
    }
  } else {
    # create output dir for this group
    output_dir <- file.path(output_dir, "sample", germline_dir, species, "vdj")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # copy the files from the original references dir
    file.copy(list.files(references_dir, full.names = TRUE), output_dir, recursive = TRUE)

    # get the loci files
    loci_files <- list.files(output_dir, paste0(loci, collapse = "|"), full.names = TRUE)

    # for each locus, and segment, build the reference and save it
    for (locus in loci) {
      locus_file <- loci_files[grep(locus, loci_files)]
      if (length(locus_file) == 0) {
        next
      }
      # get the locus genotype
      genotypes_locus <- genotypes[grepl(locus, genotypes$gene), ]

      for (seg in c("V", "D", "J")) {
        call <- paste0(locus, seg)
        # check that there are genes for the locus and segment
        if (any(grepl(call, genotypes_locus$gene))) {
          seg_references <- references[[locus]][[seg]]
          # get the genotype reference
          seg_ref_group <- .genotyped_reference(genotypes_locus[grepl(call, genotypes_locus$gene), ], seg_references)
          # write the genotype reference
          seg_ref_file <- locus_file[grepl(paste0(locus, seg, "\\."), basename(locus_file))]
          if (length(seg_ref_file) != 1) {
            stop(sprintf("Expected exactly one reference file for %s, found %d.", call, length(seg_ref_file)))
          }
          tigger::writeFasta(seg_ref_group, file.path(output_dir, basename(seg_ref_file)))
          if (!quiet) {
            message(sprintf("%s Genotyped Reference database: [%s](%s)\n", call, seg_ref_file, seg_ref_file))
          }
        } else {
          next
        }
      }
    }
  }
}

#' Plot genotype
#'
#' \code{plot_genotype} plots a genotype.
#'
#' @param    genotype       data frame of genotype information.
#' @param    allele_column  name of the column containing alleles.
#' @param    facet_by       column name to facet the plot by.
#' @param    gene_sort      method to sort genes ("name" or "position").
#' @param    text_size      size of the text in the plot.
#' @param    silent         if \code{TRUE}, do not plot the result.
#' @param    ...            additional arguments to pass to \code{theme()}.
#'
#' @return   A \code{ggplot} object.
#'
#' @export
plot_genotype <- function(genotype, allele_column = "alleles", facet_by = NULL, gene_sort = c("name", "position"),
                          text_size = 12, silent = FALSE, ...) {
  # Check arguments
  gene_sort <- match.arg(gene_sort)

  # Split genes' alleles into their own rows
  alleles <- strsplit(genotype[[allele_column]], ",")
  geno2 <- genotype[rep(seq_len(nrow(genotype)), lengths(alleles)), , drop = FALSE]
  geno2$alleles <- trimws(unlist(alleles, use.names = FALSE))

  # Set the gene order
  geno2$gene <- factor(geno2$gene,
    levels = rev(sortAlleles(unique(geno2$gene), method = gene_sort))
  )

  # Create the base plot
  p <- ggplot(geno2, aes(
    x = !!rlang::sym("gene"),
    fill = !!rlang::sym("alleles")
  )) +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = text_size),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    geom_bar(position = "fill") +
    coord_flip() +
    xlab("Gene") +
    ylab("") +
    scale_fill_hue(name = "Allele", h = c(0, 270), h.start = 10)

  # Plot, with facets by SUBJECT if that column is present
  if (!is.null(facet_by)) {
    p <- p + facet_grid(paste0(".~", facet_by))
  }

  # Add additional theme elements
  p <- p + do.call(theme, list(...))

  # Plot
  if (!silent) {
    plot(p)
  }

  invisible(p)
}
