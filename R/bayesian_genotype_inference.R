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
      write.table(genotypes_genotypeby, file = file.path(out_dir, paste0(val, "_", out_fn)), row.names = FALSE, quote = FALSE)
    }
  } else {
    write.table(genotypes, file = file.path(out_dir, out_fn), row.names = FALSE, quote = FALSE)
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
    alleles <- if (genotype[i, "genotyped_alleles"] == "") genotype[i, "alleles"] else genotype[i, "genotyped_alleles"]
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
      output_dir_genotypeby <- file.path(output_dir, as.character(val), species, "vdj")
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
            seg_ref_group <- .genotyped_reference(genotypes_genotypeby_locus[grepl(call, genotypes_genotypeby$gene), ], seg_references)
            # write the genotype reference
            seg_ref_file <- grep(seg, locus_file, value = TRUE)
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
    output_dir <- file.path(output_dir, species, "vdj")
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
          seg_ref_file <- grep(seg, locus_file, value = TRUE)
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
