#' Create a genotype project template
#'
#' Create the directory skeleton for a genotype project from the package's
#' bundled RStudio project template. The template serves all genotype methods
#' (see [genotype_methods()]); the method is selected at render time through the
#' report's `method` parameter.
#'
#' @param path Path to the directory where the project will be created.
#' @param ... Additional arguments (currently unused).
#' @export
genotype_project <- function(path, ...) {
  skeleton_dir <- file.path(
    system.file(package = "enchantr"), "rstudio",
    "templates", "project",
    "genotype_project_files"
  )
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Genotype method registry
# ---------------------------------------------------------------------------
#' Genotype inference for one repertoire using TIgGER's Bayesian method
#' @param d Repertoire records for a single group.
#' @param ctx Inference context built by `infer_genotype()`.
#' @return A TIgGER genotype data.frame.
#' @keywords internal
.infer_bayesian <- function(d, ctx) {
  tigger::inferGenotypeBayesian(
    d,
    find_unmutated = ctx$find_unmutated,
    germline_db = ctx$seg_references,
    v_call = ctx$call_col,
    seq = ctx$seq,
    priors = ctx$priors,
    genotyped_alleles = TRUE
  )
}

#' Genotype inference for one repertoire using TIgGER's fraction method
#' @inheritParams .infer_bayesian
#' @return A TIgGER genotype data.frame.
#' @keywords internal
.infer_fraction <- function(d, ctx) {
  tigger::inferGenotype(
    d,
    find_unmutated = ctx$find_unmutated,
    germline_db = ctx$seg_references,
    v_call = ctx$call_col,
    seq = ctx$seq,
    fraction_to_explain = ctx$fraction_to_explain,
    gene_cutoff = ctx$gene_cutoff
  )
}

#' PIgLET's bundled allele threshold table
#'
#' `piglet::inferGenotypeAllele()` loads its default with a bare
#' `data(allele_threshold_table)`, which only resolves when piglet is attached;
#' enchantr calls piglet with `::`, so load it explicitly.
#'
#' @return PIgLET's default allele threshold table.
#' @keywords internal
.piglet_default_thresholds <- function() {
  e <- new.env(parent = emptyenv())
  utils::data("allele_threshold_table", package = "piglet", envir = e)
  e$allele_threshold_table
}

#' Genotype inference for one repertoire using PIgLET's allele-based method
#' @inheritParams .infer_bayesian
#' @return A PIgLET genotype data.frame (long: one row per allele).
#' @keywords internal
.infer_allele_based <- function(d, ctx) {
  thresholds <- ctx$allele_threshold_table
  if (is.null(thresholds)) {
    thresholds <- .piglet_default_thresholds()
  }
  piglet::inferGenotypeAllele(
    d,
    allele_threshold_table = thresholds,
    call = ctx$call_col,
    seq = ctx$seq,
    find_unmutated = ctx$find_unmutated,
    germline_db = ctx$seg_references,
    default_allele_threshold = ctx$default_threshold,
    asc_annotation = ctx$asc_annotation,
    translate_to_asc = ctx$translate_to_asc,
    depth_adjusted_threshold = ctx$depth_adjusted_threshold,
    z_score_threshold = ctx$z_score_threshold,
    single_assignment = FALSE
  )
}

#' Strip the gene prefix from allele calls
#'
#' `IGHV1-2*02` -> `02`. Calls that carry no `*` are returned unchanged.
#'
#' @param x Character vector of allele calls.
#' @return Character vector of bare allele numbers.
#' @keywords internal
.bare_allele <- function(x) sub(".*\\*", "", x)

#' Normalize TIgGER genotype output onto the shared genotype schema
#'
#' `inferGenotypeBayesian()` scores its call with `k_diff`, normalized to
#' `confidence`. `inferGenotype()` reports the genotype in `alleles` and has no
#' confidence measure, so it emits no `confidence` column at all -- an all-NA
#' column would read as a failed computation rather than an absent metric.
#'
#' @param g Genotype data.frame as returned by TIgGER.
#' @param method Method name.
#' @return A data.frame following the shared genotype schema.
#' @keywords internal
.normalize_tigger <- function(g, method) {
  if (is.null(g) || nrow(g) == 0) {
    return(NULL)
  }

  g$candidate_alleles <- g$alleles
  if (identical(method, "fraction")) {
    # no separate genotyped column: `alleles` is the genotype
    g$genotyped_alleles <- g$alleles
  } else {
    # TIgGER returns k_diff as character; the confidence plot needs a numeric scale
    g$confidence <- as.numeric(g$k_diff)
    g$k_diff <- NULL
  }
  g$alleles <- NULL

  if (!"note" %in% colnames(g)) {
    g$note <- ""
  }
  g
}

#' Normalize PIgLET genotype output onto the shared genotype schema
#'
#' Pivots PIgLET's long output (one row per candidate allele, scored by
#' `z_score`) to one row per gene. Genotype membership is PIgLET's own
#' `in_genotype` column, which already applies `z_score_threshold` and any depth
#' adjustment -- re-deriving it here would silently discard both settings.
#' Rows whose `allele` is itself a comma-separated
#' call are assignment bookkeeping and are dropped, as are genes with no allele
#' called. `depth` is PIgLET's repertoire-wide depth, so `total` is recomputed per
#' gene to match TIgGER. `confidence` carries each genotyped allele's own
#' `z_score`, comma-separated and aligned with `genotyped_alleles` -- PIgLET
#' scores alleles, not genes, so no gene-level summary is invented here.
#'
#' @param g Genotype data.frame as returned by PIgLET.
#' @param method Method name (unused; kept for a uniform normalize signature).
#' @return A data.frame following the shared genotype schema.
#' @keywords internal
.normalize_piglet <- function(g, method) {
  if (is.null(g) || nrow(g) == 0) {
    return(NULL)
  }

  g <- as.data.frame(g, stringsAsFactors = FALSE)
  g <- g[!grepl(",", g$allele), , drop = FALSE]
  if (nrow(g) == 0) {
    return(NULL)
  }
  g$.in_genotype <- !is.na(g$in_genotype) & g$in_genotype

  collapse <- function(x) paste(x, collapse = ",")

  out <- g %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      genotyped_alleles = collapse(.bare_allele(allele[.in_genotype])),
      candidate_alleles = collapse(.bare_allele(allele)),
      # counts align with candidate_alleles, matching TIgGER's convention.
      counts = collapse(round(count, 2)),
      total = round(sum(count), 2),
      depth = dplyr::first(depth),
      confidence = collapse(round(z_score[.in_genotype], 2)),
      note = "",
      .groups = "drop"
    ) %>%
    dplyr::filter(nzchar(genotyped_alleles)) %>%
    as.data.frame(stringsAsFactors = FALSE)

  if (nrow(out) == 0) NULL else out
}

#' Plot a genotype without a confidence scale
#' @param g Normalized genotype data.frame.
#' @param ... Passed to `tigger::plotGenotype()`.
#' @return A ggplot object.
#' @keywords internal
.plot_plain <- function(g, ...) {
  tigger::plotGenotype(g, allele_col = "genotyped_alleles", silent = TRUE, ...)
}

#' Plot a genotype shaded by the method's confidence measure
#' @param g Normalized genotype data.frame.
#' @param ... Passed to `tigger::plotGenotypeConfidence()`.
#' @return A ggplot object.
#' @keywords internal
.plot_confidence <- function(g, ...) {
  tigger::plotGenotypeConfidence(
    g,
    confidence_col = "confidence",
    allele_col = "genotyped_alleles",
    silent = TRUE,
    ...
  )
}

# Per-method functions, parameters, extra columns and report text.
.GENOTYPE_METHODS <- list(
  bayesian = list(
    label = "TIgGER Bayesian",
    package = "tigger",
    params = c("v_priors", "d_priors", "j_priors"),
    extra_cols = c("kh", "kd", "kt", "kq"),
    infer = .infer_bayesian,
    normalize = .normalize_tigger,
    plot = .plot_confidence,
    description = paste(
      "The `inferGenotypeBayesian` function applies a Bayesian framework to infer subject-specific germline genotypes from observed V(D)J gene assignments.",
      "Allele usage is modeled as a multinomial distribution with a Dirichlet prior, empirically parameterized from high-coverage reference genotypes.",
      "",
      "Model certainty is summarized by `k_diff`, the log10 ratio between the highest- and second-highest-likelihood zygosity models. Larger `k_diff` values indicate clearer separation between competing models and therefore stronger support for the inferred genotype. Lower `k_diff` values, especially values closer to zero, indicate that the leading models are poorly separated. In practice, low `k_diff` often reflects a low number of informative sequences for that gene, so confidence in the call is lower.",
      "",
      "This probabilistic approach incorporates prior expectations of allele frequencies within each gene family and accounts for sequencing depth variation, ambiguous V assignments, and low-frequency alleles.",
      sep = "\n"
    ),
    columns = paste(
      "- `confidence`: `k_diff`, the log10 ratio between the best- and second-best-supported zygosity models.",
      "- `kh`, `kd`, `kt`, `kq`: log10 likelihoods for the homozygous, heterozygous, trizygous, and quadrozygous genotype models.",
      "",
      "Larger `confidence` values indicate stronger separation between competing genotype models; values closer to zero indicate weaker support for the selected genotype.",
      sep = "\n"
    )
  ),
  fraction = list(
    label = "TIgGER fraction",
    package = "tigger",
    params = c("fraction_to_explain", "gene_cutoff"),
    extra_cols = character(0),
    infer = .infer_fraction,
    normalize = .normalize_tigger,
    plot = .plot_plain,
    description = paste(
      "The `inferGenotype` function infers the genotype using TIgGER's cumulative-fraction method.",
      "For each gene, alleles are ranked by the number of sequences assigned to them, and alleles are added to the genotype until they cumulatively explain `fraction_to_explain` of the sequences for that gene. Alleles falling below `gene_cutoff` are discarded as noise.",
      "",
      "This method is deterministic and easy to reason about, but unlike the Bayesian method it produces **no per-gene confidence measure** -- the genotype plot below is therefore unshaded. Genes supported by few sequences are not distinguished from well-supported ones, so the `sequence_count` and `clone_count` columns should be read carefully.",
      sep = "\n"
    ),
    columns = paste(
      "The fraction method does not score competing genotype models, so it has no confidence measure and the genotype table carries **no `confidence` column**. Judge support from `sequence_count` and `clone_count` instead.",
      sep = "\n"
    )
  ),
  allele_based = list(
    label = "PIgLET allele based",
    package = "piglet",
    params = c("allele_thresholds_db", "default_threshold", "asc_annotation", "translate_to_asc", "z_score_threshold"),
    extra_cols = "depth",
    infer = .infer_allele_based,
    normalize = .normalize_piglet,
    # `confidence` is per-allele here, so it cannot feed .plot_confidence,
    # which expects one numeric score per gene.
    plot = .plot_plain,
    description = paste(
      "The `inferGenotypeAllele` function from PIgLET infers the genotype using allele-specific thresholds.",
      "Rather than applying one rule to every allele, each allele carries its own population-derived presence threshold. From the allele's count, the repertoire depth, and that threshold, a standardized score (`z_score`) is computed for the allele being present in the genotype. Alleles whose `z_score` reaches `z_score_threshold` are called into the genotype; raising it above its default of 0 demands more evidence per allele and yields a smaller, more conservative genotype.",
      "",
      "Because thresholds are set per allele, this method handles alleles that are systematically hard to detect -- for example those that are frequently mis-assigned or that sit in low-coverage regions -- more even-handedly than a single global cutoff. Alleles absent from the threshold table fall back to `default_threshold`.",
      "",
      "Each allele's threshold is first raised to at least `max(default_threshold, 1 / depth)`, so no allele is expected at a frequency below a single observed sequence. This matters most in shallow repertoires, where a population-derived threshold can otherwise sit below the resolution the data actually supports; in deep repertoires `1 / depth` is small and the adjustment rarely changes the call. The report always applies it, so it is not a parameter.",
      sep = "\n"
    ),
    columns = paste(
      "- `confidence`: the `z_score` of each genotyped allele, comma-separated and in the same order as `genotyped_alleles`. PIgLET scores alleles individually, so there is one value per allele rather than one per gene.",
      "- `depth`: the repertoire depth PIgLET used as the denominator of every `z_score`. It is the sum of the ambiguity-weighted allele `counts` across the whole locus, not a count of sequences: an ambiguous call is split across the alleles it names and down-weighted, and alleles with no observation at all receive a small pseudo-count. `depth` is therefore fractional, is usually smaller than the number of sequences analysed, and is identical on every row of a locus.",
      "",
      "Alleles were called into the genotype when their individual `z_score` reached `z_score_threshold`. See the input parameter table above for the value used.",
      sep = "\n"
    )
  )
)

#' Available genotype inference methods
#'
#' @return Character vector of supported method names.
#' @export
genotype_methods <- function() names(.GENOTYPE_METHODS)

#' Human-readable name of a genotype method
#'
#' Used for the report subtitle and the `METHOD>` line of the command log, so the
#' label for each method is defined in one place.
#'
#' @param method The method name.
#' @return A character scalar, e.g. "TIgGER Bayesian".
#' @export
genotype_method_label <- function(method) .genotype_method(method)$label

#' Look up a genotype method, with a helpful error for unknown names
#' @param method Method name.
#' @return The method's registry entry.
#' @keywords internal
.genotype_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || !method %in% names(.GENOTYPE_METHODS)) {
    stop(
      "Unknown genotype method: ", paste(method, collapse = ", "),
      ". Available methods: ", paste(genotype_methods(), collapse = ", "), "."
    )
  }
  .GENOTYPE_METHODS[[method]]
}

#' Validate report parameters against the selected genotype method
#'
#' The report template must declare the parameters of every method, because R
#' Markdown's `params` block is static -- it cannot depend on the value of
#' `method`. This drops the parameters belonging to the other methods at run
#' time, so the rendered report only advertises the ones that actually applied.
#'
#' @param params The report's `params` list.
#' @param method The selected method name.
#' @return `params` with the other methods' parameters removed.
#' @export
check_genotype_params <- function(params, method) {
  spec <- .genotype_method(method)
  foreign <- setdiff(
    unlist(lapply(.GENOTYPE_METHODS, `[[`, "params"), use.names = FALSE),
    spec$params
  )
  ignored <- intersect(names(params), foreign)
  if (length(ignored) > 0) {
    message(
      "Genotype method '", method, "' does not use these parameters; ignoring: ",
      paste(ignored, collapse = ", "), "."
    )
  }

  params[setdiff(names(params), ignored)]
}

#' Was an optional file parameter actually supplied?
#'
#' Report parameters arrive as NULL, `NA`, `""` or the literal string `"NULL"`
#' when left unset, depending on how the report was invoked.
#'
#' @param path A report parameter naming a file.
#' @return TRUE when `path` is a usable path.
#' @export
has_path <- function(path) {
  !is.null(path) && length(path) == 1 && !is.na(path) &&
    nzchar(path) && !identical(path, "NULL")
}

#' Read an allele threshold table
#'
#' Loads the allele-specific presence thresholds used by the `"allele_based"`
#' method. When no path is given, NULL is returned so that
#' `piglet::inferGenotypeAllele()` falls back to its own bundled thresholds.
#'
#' @param path Path to a TSV of allele thresholds, or NULL/"NULL"/"" for none.
#' @return A data.frame of thresholds, or NULL.
#' @export
read_allele_thresholds <- function(path) {
  if (!has_path(path)) {
    return(NULL)
  }
  utils::read.delim(path, stringsAsFactors = FALSE)
}

#' Explanatory text for a genotype method
#'
#' The report's method-varying prose lives in `.GENOTYPE_METHODS` rather than in
#' the template, so that a single template can serve every method.
#'
#' @param method The selected method name.
#' @param section Which block of text to return: `"description"` explains how the
#'   method infers the genotype; `"columns"` documents the method-specific columns
#'   of the genotype table.
#' @return A markdown string.
#' @export
genotype_method_text <- function(method, section = c("description", "columns")) {
  .genotype_method(method)[[match.arg(section)]]
}

#' Plot an inferred genotype
#'
#' Renders the genotype table produced by [infer_genotype()], using the plot that
#' suits the method: methods that report a confidence measure are shaded by it,
#' methods that do not are drawn as a plain presence/absence map.
#'
#' @param genotypes A normalized genotype data.frame (see [infer_genotype()]).
#' @param method The genotype method that produced `genotypes`.
#' @param ... Passed through to the underlying TIgGER plotting function.
#' @return A ggplot object.
#' @export
plot_genotype <- function(genotypes, method = "bayesian", ...) {
  spec <- .genotype_method(method)
  spec$plot(genotypes, ...)
}

## wrapper to infer genotype by segment, and by genotypeby column and return the results
#' Infer a genotype for a segment
#'
#' Runs genotype inference for a single segment (V/D/J) using the requested
#' `method`, and normalizes the result onto a shared schema so that every method
#' can feed the same downstream table, plot, and reference-generation code.
#'
#' Regardless of method, the returned data.frame has one row per gene and
#' contains at least `gene`, `genotyped_alleles`, `candidate_alleles`, `counts`,
#' `total`, `confidence`, and `note`. `confidence` carries whatever the method
#' uses to score its call (`k_diff` for `"bayesian"`, `z_score` for
#' `"allele_based"`, `NA` for `"fraction"`, which has no confidence measure).
#'
#' Optionally the inference can be run per-group by providing `genotypeby`, in
#' which case results are tagged with the group value and row-bound.
#'
#' @param db A data.frame or tibble containing rearrangement records (AIRR format).
#' @param seg Character scalar. Segment short name, e.g. "v", "d", or "j".
#' @param loci Character vector. Loci to extract reference alleles from when
#'   `references` is provided (e.g. c("IGH")).
#' @param method Genotype inference method; one of [genotype_methods()].
#'   `"bayesian"` uses `tigger::inferGenotypeBayesian()`, `"fraction"` uses
#'   `tigger::inferGenotype()`, and `"allele_based"` uses
#'   `piglet::inferGenotypeAllele()`.
#' @param genotypeby Optional character column name in `db` to run inference per-group.
#' @param references Optional nested list of reference sequences (as produced by dowser/readIMGT). If provided,
#'   alleles for the selected `seg` are collected from each `locus` in `loci` and
#'   passed to the method via `germline_db`.
#' @param find_unmutated Logical. Use `references` to restrict inference to unmutated sequences.
#' @param single_assignments Logical. If TRUE, only sequences with single allele
#'   assignments (no commas) are used for inference. Applied once, up front, so
#'   that every method sees the same input.
#' @param seq Column name holding alignment sequences (default: "sequence_alignment").
#' @param priors `"bayesian"` only. Numeric vector of prior probabilities for the
#'   genotype model; see `tigger::inferGenotypeBayesian()`.
#' @param fraction_to_explain `"fraction"` only. Cumulative fraction of a gene's
#'   sequences the genotyped alleles must explain.
#' @param gene_cutoff `"fraction"` only. Alleles below this frequency are discarded.
#' @param allele_threshold_table `"allele_based"` only. Data.frame of alleles and
#'   their population-derived presence thresholds.
#' @param default_threshold `"allele_based"` only. Threshold used for alleles absent
#'   from `allele_threshold_table`.
#' @param asc_annotation `"allele_based"` only. Are the allele calls annotated with
#'   allele similarity clusters?
#' @param translate_to_asc `"allele_based"` only. Collapse identical alleles for inference.
#' @param depth_adjusted_threshold `"allele_based"` only. Scale allele thresholds by
#'   the repertoire depth.
#' @param z_score_threshold `"allele_based"` only. Minimum `z_score` for an allele to
#'   be called into the genotype.
#' @return A normalized genotype data.frame (see Details), or NULL if inference
#'   produced no result.
#' @export
infer_genotype <- function(
  db,
  seg,
  loci,
  method = "bayesian",
  genotypeby = NULL,
  references = NULL,
  single_assignments = FALSE,
  find_unmutated = FALSE,
  seq = "sequence_alignment",
  priors = c(0.6, 0.4, 0.4, 0.35, 0.25, 0.25, 0.25, 0.25, 0.25),
  fraction_to_explain = 0.875,
  gene_cutoff = 1e-04,
  allele_threshold_table = NULL,
  default_threshold = 1e-04,
  asc_annotation = FALSE,
  translate_to_asc = FALSE,
  depth_adjusted_threshold = TRUE,
  z_score_threshold = 0
) {
  spec <- .genotype_method(method)
  if (!requireNamespace(spec$package, quietly = TRUE)) {
    stop(
      "Genotype method '", method, "' requires the '", spec$package,
      "' package, which is not installed."
    )
  }

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

  # Everything the method's inference function needs, assembled once.
  ctx <- list(
    call_col = call_col,
    seq = seq,
    find_unmutated = find_unmutated,
    seg_references = seg_references,
    priors = priors,
    fraction_to_explain = fraction_to_explain,
    gene_cutoff = gene_cutoff,
    allele_threshold_table = allele_threshold_table,
    default_threshold = default_threshold,
    asc_annotation = asc_annotation,
    translate_to_asc = translate_to_asc,
    depth_adjusted_threshold = depth_adjusted_threshold,
    z_score_threshold = z_score_threshold
  )

  # A failed group yields NULL with a warning instead of stopping the run.
  run_one <- function(d) {
    g <- try(spec$infer(d, ctx), silent = TRUE)
    if (inherits(g, "try-error")) {
      warning(
        "Genotype inference failed for segment '", seg, "' using method '", method,
        "': ", conditionMessage(attr(g, "condition")),
        call. = FALSE
      )
      return(NULL)
    }
    spec$normalize(g, method)
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

#' Record one filtering step
#'
#' Every filter reports rows in and rows out, so a step that removes sequences is
#' visible in the report rather than only in the final count.
#'
#' @param steps A data.frame of previous steps, or NULL.
#' @param label What the step did.
#' @param rows_in,rows_out Row counts before and after.
#' @return The steps table with one row appended.
#' @export
count_rows_step <- function(steps, label, rows_in, rows_out) {
  rbind(steps, data.frame(
    step = label, rows_in = rows_in, rows_out = rows_out,
    dropped = rows_in - rows_out, stringsAsFactors = FALSE
  ))
}

#' The column a threshold table is looked up by
#'
#' PIgLET matches thresholds on `asc_allele` when the calls are ASC names and on
#' `allele` otherwise. The table carries both, so it is also the crosswalk
#' between the two naming schemes -- which is why the naming is never guessed
#' from the data.
#'
#' @param thresholds A threshold table.
#' @param asc_annotation Logical. TRUE when the calls in the repertoire are ASC names.
#' @return The name of the column to match allele calls against.
#' @keywords internal
.threshold_key_column <- function(thresholds, asc_annotation = FALSE) {
  key <- if (isTRUE(asc_annotation)) "asc_allele" else "allele"
  if (!key %in% colnames(thresholds)) {
    stop(
      "The allele threshold table has no '", key, "' column, so allele calls cannot be ",
      "looked up in it. Columns present: ", paste(colnames(thresholds), collapse = ", "), ". ",
      if (identical(key, "asc_allele")) {
        "With asc_annotation = TRUE the table must carry the ASC names."
      } else {
        "With asc_annotation = FALSE the table must carry the IMGT-style allele names."
      }
    )
  }
  key
}

#' How much of the repertoire the threshold table actually covers
#'
#' Reports, per allele observed in the repertoire, whether a threshold was found
#' for it. An allele with no entry silently inherits `default_threshold`, which
#' discards the population-derived value the method is supposed to use -- so the
#' coverage is reported, and zero coverage is an error rather than a fallback.
#'
#' @param db A data.frame of repertoire records.
#' @param call_col The allele call column, e.g. `v_call`.
#' @param thresholds A threshold table, or NULL.
#' @param asc_annotation Logical, see [.threshold_key_column()].
#' @param single_assignments Logical. The rule the inference applied; when FALSE an
#'   ambiguous call contributes every allele it lists.
#' @return A data.frame with one row per observed allele: `allele`,
#'   `sequence_count`, `matched`, and the `threshold` when one was found.
#' @export
allele_threshold_coverage <- function(db, call_col, thresholds,
                                      asc_annotation = FALSE, single_assignments = FALSE) {
  if (is.null(thresholds) || nrow(thresholds) == 0 || is.null(db) || nrow(db) == 0 ||
    !call_col %in% colnames(db)) {
    return(data.frame())
  }
  key <- .threshold_key_column(thresholds, asc_annotation)

  calls <- db[[call_col]]
  calls <- calls[!is.na(calls) & nzchar(calls)]
  if (single_assignments) {
    calls <- calls[!grepl(",", calls)]
  } else {
    calls <- trimws(unlist(strsplit(calls, ",")))
  }
  if (length(calls) == 0) {
    return(data.frame())
  }

  observed <- as.data.frame(table(allele = calls), stringsAsFactors = FALSE)
  colnames(observed) <- c("allele", "sequence_count")
  lookup <- thresholds[match(observed$allele, thresholds[[key]]), , drop = FALSE]
  observed$matched <- !is.na(lookup[[key]])
  observed$threshold <- if ("threshold" %in% colnames(thresholds)) lookup$threshold else NA
  observed[order(-observed$sequence_count), , drop = FALSE]
}

#' Summarize threshold coverage into one row
#'
#' @param coverage Output of [allele_threshold_coverage()].
#' @return A one-row data.frame: alleles observed, matched, defaulted, and the
#'   matched fraction.
#' @export
summarize_threshold_coverage <- function(coverage) {
  if (is.null(coverage) || nrow(coverage) == 0) {
    return(data.frame(
      alleles_observed = 0L, alleles_matched = 0L, alleles_defaulted = 0L,
      sequences_matched = 0L, sequences_defaulted = 0L, coverage = NA_real_
    ))
  }
  data.frame(
    alleles_observed = nrow(coverage),
    alleles_matched = sum(coverage$matched),
    alleles_defaulted = sum(!coverage$matched),
    sequences_matched = sum(coverage$sequence_count[coverage$matched]),
    sequences_defaulted = sum(coverage$sequence_count[!coverage$matched]),
    coverage = sum(coverage$matched) / nrow(coverage)
  )
}

#' Provenance of the allele threshold table in use
#'
#' A genotype is only reproducible if the thresholds behind it can be identified.
#' PIgLET's bundled table carries no release attribute, so the checksum of the
#' table actually used is what pins the run.
#'
#' @param thresholds The threshold table in use.
#' @param path Path it was read from, or NULL when PIgLET's bundled table is used.
#' @param release Release label supplied by the caller, or NULL when unspecified.
#' @return A one-row data.frame of provenance fields.
#' @export
allele_threshold_provenance <- function(thresholds, path = NULL, release = NULL) {
  checksum <- NA_character_
  if (!is.null(thresholds) && nrow(thresholds) > 0) {
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    utils::write.table(thresholds, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
    checksum <- unname(tools::md5sum(tmp))
  }
  data.frame(
    source = if (has_path(path)) path else "piglet bundled allele_threshold_table",
    release = if (!is.null(release) && nzchar(release)) release else "unspecified",
    n_alleles = if (is.null(thresholds)) 0L else nrow(thresholds),
    columns = if (is.null(thresholds)) "" else paste(colnames(thresholds), collapse = ","),
    md5 = checksum,
    piglet_version = tryCatch(as.character(utils::packageVersion("piglet")), error = function(e) NA_character_),
    recorded = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
}

#' Summarize support for genotype genes
#'
#' Internal helper that counts productive sequences and unique clones supporting
#' each genotyped allele. The counts follow the same assignment rule the
#' inference used, so that the reported support and the genotype are computed
#' over the same sequences: with `single_assignments = TRUE` only unambiguous
#' calls are counted, otherwise a sequence supports every allele its call lists
#' and so contributes to more than one allele.
#'
#' @param db A data.frame/tibble of productive repertoire records.
#' @param genotypeby Optional grouping column name used for per-group genotype
#'   inference.
#' @param clone_id_column Optional clone identifier column name.
#' @param single_assignments Logical. The rule `infer_genotype()` applied to the
#'   same repertoire; pass the value used there, not a new choice.
#' @return A data.frame with `gene`, `allele`, `sequence_count`,
#'   `clone_count`, and the grouping column when `genotypeby` is provided.
#' @keywords internal
.genotype_support_counts <- function(db, genotypeby = NULL, clone_id_column = NULL,
                                    single_assignments = FALSE) {
  if (is.null(db) || nrow(db) == 0) {
    return(data.frame())
  }

  count_by_segment <- lapply(c("v", "d", "j"), function(seg) {
    call_col <- paste0(seg, "_call")
    if (!call_col %in% colnames(db)) {
      return(NULL)
    }

    seg_db <- db[!is.na(db[[call_col]]) & nzchar(db[[call_col]]), , drop = FALSE]
    if (single_assignments) {
      # the inference saw only unambiguous calls, so count only those
      seg_db <- seg_db[!grepl(",", seg_db[[call_col]]), , drop = FALSE]
    } else {
      # an ambiguous call supports every allele it lists
      calls <- strsplit(seg_db[[call_col]], ",")
      seg_db <- seg_db[rep(seq_len(nrow(seg_db)), lengths(calls)), , drop = FALSE]
      seg_db[[call_col]] <- trimws(unlist(calls))
    }
    if (nrow(seg_db) == 0) {
      return(NULL)
    }

    seg_db$gene <- alakazam::getGene(seg_db[[call_col]], strip_d = FALSE)
    seg_db$allele <- .bare_allele(seg_db[[call_col]])
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
#' @param method The genotype method that produced `genotypes`. Selects which
#'   method-specific columns are placed in the presentation order; it is not
#'   written as a column.
#' @return A data.frame with added `sequence_count`/`clone_count` and reordered
#'   columns.
#' @export
augment_genotype_table <- function(genotypes, db_support, genotypeby = NULL, clone_id_column = NULL,
                                  single_assignments = FALSE, method) {
  spec <- .genotype_method(method)
  if (is.null(genotypes) || nrow(genotypes) == 0) {
    return(genotypes)
  }

  support_counts <- .genotype_support_counts(
    db = db_support,
    genotypeby = genotypeby,
    clone_id_column = clone_id_column,
    single_assignments = single_assignments
  )
  group_cols <- c(if (!is.null(genotypeby) && genotypeby %in% colnames(genotypes)) genotypeby, "gene")
  has_clone_counts <- nrow(support_counts) > 0 && !all(is.na(support_counts$clone_count))

  # one row per (genotype row, genotyped allele), joined to its support counts
  alleles <- strsplit(genotypes$genotyped_alleles, ",")
  alleles[lengths(alleles) == 0] <- ""
  per_allele <- data.frame(
    .row = rep(seq_len(nrow(genotypes)), lengths(alleles)),
    allele = trimws(unlist(alleles)),
    stringsAsFactors = FALSE
  )
  for (col in group_cols) {
    per_allele[[col]] <- as.character(genotypes[[col]])[per_allele$.row]
    if (nrow(support_counts) > 0) support_counts[[col]] <- as.character(support_counts[[col]])
  }
  if (nrow(support_counts) > 0) {
    per_allele <- dplyr::left_join(per_allele, support_counts, by = c(group_cols, "allele"))
  } else {
    per_allele$sequence_count <- NA_integer_
    per_allele$clone_count <- NA_integer_
  }

  # an allele the inference called but the counts never saw contributes 0
  collapse_counts <- function(x) {
    out <- tapply(x, per_allele$.row, function(v) paste(ifelse(is.na(v), "0", as.character(v)), collapse = ","))
    as.character(out[as.character(seq_len(nrow(genotypes)))])  # tapply returns a 1-d array
  }
  genotypes$sequence_count <- collapse_counts(per_allele$sequence_count)
  genotypes$clone_count <- if (has_clone_counts) collapse_counts(per_allele$clone_count) else NA_character_
  no_alleles <- is.na(genotypes$genotyped_alleles) | !nzchar(genotypes$genotyped_alleles)
  genotypes$sequence_count[no_alleles] <- NA_character_
  genotypes$clone_count[no_alleles] <- NA_character_

  # shared schema first, then the columns only this method produces
  core_cols <- c(
    if (!is.null(genotypeby) && genotypeby %in% colnames(genotypes)) genotypeby,
    "gene",
    "genotyped_alleles",
    "sequence_count",
    "clone_count",
    "confidence",
    "note",
    "candidate_alleles",
    "counts",
    "total"
  )
  ordered_cols <- intersect(c(core_cols, spec$extra_cols), colnames(genotypes))
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
  write_one <- function(df, fn) {
    write.table(df, file = file.path(out_dir, fn), row.names = FALSE, quote = FALSE, sep = "\t")
  }
  if (is.null(genotypeby)) {
    write_one(genotypes, out_fn)
  } else {
    for (val in unique(genotypes[[genotypeby]])) {
      write_one(genotypes[genotypes[[genotypeby]] == val, ], paste0(val, "_", out_fn))
    }
  }
}

#' Generate genotype-informed reference FASTA files
#'
#' For a combined genotype table (optionally containing a grouping column
#' given by `genotypeby`), build per-locus, per-segment genotype-aware
#' reference FASTA files via \code{tigger::genotypeFasta(..., include_unseen = TRUE)}
#' (genes absent from the genotype keep all their reference alleles; genes
#' present are restricted to their genotyped alleles). When `genotypeby` is
#' provided a separate output tree is created per group, otherwise a single
#' "sample" tree is written.
#'
#' The baseline files from `references_dir` are copied into the target location
#' and the per-locus FASTA files are overwritten with the genotyped sequences.
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
  # tigger::genotypeFasta expects `alleles`
  if (!"alleles" %in% colnames(genotypes) && "candidate_alleles" %in% colnames(genotypes)) {
    genotypes$alleles <- genotypes$candidate_alleles
  }

  # Treat the non-grouped case as a single "sample" group so the loop is shared.
  groups <- if (!is.null(genotypeby)) unique(genotypes[[genotypeby]]) else "sample"

  for (val in groups) {
    group_genotypes <- if (is.null(genotypeby)) genotypes else genotypes[genotypes[[genotypeby]] == val, ]

    out_dir <- file.path(output_dir, as.character(val), germline_dir, species, "vdj")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(references_dir, full.names = TRUE), out_dir, recursive = TRUE)
    loci_files <- list.files(out_dir, paste0(loci, collapse = "|"), full.names = TRUE)

    for (locus in loci) {
      locus_file <- loci_files[grep(locus, loci_files)]
      if (length(locus_file) == 0) {
        next
      }
      group_locus <- group_genotypes[grepl(locus, group_genotypes$gene), ]

      for (seg in c("V", "D", "J")) {
        call <- paste0(locus, seg)
        seg_references <- references[[locus]][[seg]]
        # skip segments with no genotyped genes or no reference
        if (!any(grepl(call, group_locus$gene)) || is.null(seg_references)) {
          next
        }

        seg_ref_group <- tigger::genotypeFasta(
          group_locus[grepl(call, group_locus$gene), ], seg_references,
          include_unseen = TRUE, strip_d = FALSE
        )
        seg_ref_file <- locus_file[grepl(paste0(locus, seg, "\\."), basename(locus_file))]
        if (length(seg_ref_file) != 1) {
          stop(sprintf("Expected exactly one reference file for %s, found %d.", call, length(seg_ref_file)))
        }
        tigger::writeFasta(seg_ref_group, file.path(out_dir, basename(seg_ref_file)))
        if (!quiet) {
          message(sprintf("%s Genotyped Reference database for %s: [%s](%s)\n", call, val, seg_ref_file, seg_ref_file))
        }
      }
    }
  }
}
