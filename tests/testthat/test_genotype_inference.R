IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

# The shared genotype schema, in the order augment_genotype_table() writes it.
REPORT_COLUMNS <- c(
  "gene", "genotyped_alleles", "sequence_count", "clone_count",
  "confidence", "note", "candidate_alleles", "counts", "total"
)

# The unified report, and the tigger_bayesian_genotype report kept for callers
# that still use it, which writes the Bayesian score as k_diff.
REPORTS <- list(
  genotype = list(params = list(method = "bayesian"), score = "confidence",
                  output = "sample_genotype_report.tsv"),
  tigger_bayesian_genotype = list(params = list(), score = "k_diff",
                                  output = "tigger_genotype_report.tsv")
)

test_that("every genotype method normalizes onto the shared schema", {
  skip_on_cran()
  db <- tigger::AIRRDb
  db <- db[!is.na(db$v_call), ]

  gene_sets <- list()
  for (method in genotype_methods()) {
    skip_if_not_installed(enchantr:::.genotype_method(method)$package)
    genotype <- infer_genotype(db, seg = "v", loci = "IGH", method = method, find_unmutated = FALSE)
    gene_sets[[method]] <- sort(genotype$gene)

    # `fraction` has no confidence measure and so emits no such column.
    required <- setdiff(REPORT_COLUMNS, c("sequence_count", "clone_count",
                                          if (method == "fraction") "confidence"))
    expect_true(all(required %in% names(genotype)), info = method)
    if (method == "fraction") {
      expect_null(genotype$confidence)
    } else {
      expect_false(any(is.na(genotype$confidence)), info = method)
      expect_true(all(nzchar(genotype$confidence)), info = method)
    }
    if (method == "allele_based") {
      # PIgLET scores alleles, not genes: one z_score per genotyped allele, in order.
      expect_equal(
        lengths(strsplit(genotype$confidence, ",")),
        lengths(strsplit(genotype$genotyped_alleles, ","))
      )
    }
    if (method == "bayesian") {
      expect_type(genotype$confidence, "double")
    }
    expect_false("alleles" %in% names(genotype), info = method)
    # one row per gene
    expect_false(any(duplicated(genotype$gene)), info = method)
    expect_true(all(nzchar(genotype$genotyped_alleles)), info = method)
  }

  # methods should agree on which genes are present
  for (method in names(gene_sets)) {
    expect_equal(gene_sets[[method]], gene_sets[[1]], info = method)
  }
})

test_that("unknown genotype method is rejected", {
  expect_error(infer_genotype(data.frame(), seg = "v", loci = "IGH", method = "nope"), "Unknown genotype method")
})

test_that("both Bayesian genotype reports: genotype inference 1:1", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))

  for (report in names(REPORTS)) {
    tmp_dir <- tempfile(paste0(report, "_1_1_"))
    enchantr_report(report, report_params = c(list(
      "input" = input, "imgt_db" = IMGT_URL, "species" = "human",
      "outdir" = tmp_dir, "log" = "test_allele_inference_command_log"
    ), REPORTS[[report]]$params))

    genotype <- read.delim(list.files(file.path(tmp_dir, "enchantr", "genotypes"), full.names = TRUE), sep = "\t")
    expected_cols <- sub("^confidence$", REPORTS[[report]]$score, c(REPORT_COLUMNS, "kh", "kd", "kt", "kq"))
    expect_equal(nrow(genotype), 41, info = report)
    expect_false("alleles" %in% names(genotype), info = report)
    expect_equal(names(genotype)[seq_along(expected_cols)], expected_cols, info = report)
    expect_true(all(is.na(genotype$clone_count)), info = report)
  }
})

test_that("both Bayesian genotype reports keep a novel allele", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))

  tmp_dir <- tempfile("genotype_novel_inference_")
  enchantr_report("novel_allele_inference",
    report_params = list(
      "input" = input,
      "imgt_db" = IMGT_URL,
      "species" = "human",
      "outdir" = tmp_dir,
      "pos_range" = "1:318",
      "nproc" = 1,
      "log" = "test_allele_inference_command_log"
    )
  )
  novel_report_dir <- file.path(tmp_dir, "enchantr")
  novel_db <- file.path(novel_report_dir, "db_novel")
  evidence_path <- file.path(novel_report_dir, "tigger-novel_novel_allele_evidence.rda")

  tmp_dir <- tempfile("genotype_novel_reassign_")
  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = input,
      "imgt_db" = novel_db,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "segments" = "v",
      "log" = "test_reassign_alleles_command_log"
    )
  )
  reassigned_input <- list.files(file.path(tmp_dir, "enchantr", "repertoires"), full.names = TRUE)
  expect_length(reassigned_input, 1)

  db <- read.delim(reassigned_input, sep = "\t", quote = "")
  target_call <- "IGHV1-2*02_A318G"
  target_gene <- sub("\\*.*", "", target_call)
  gene_db <- db[grepl(paste0("^", target_gene, "\\*"), db$v_call), ]
  target_db <- db[db$v_call == target_call, ]
  expect_gt(nrow(target_db), 0)

  genotype_input <- rbind(gene_db, target_db[rep(seq_len(nrow(target_db)), 20), ])
  genotype_input$sequence_id <- paste0("genotype_novel_", seq_len(nrow(genotype_input)))
  genotype_input_path <- tempfile("genotype_input_with_novel_", fileext = ".tsv")
  write.table(genotype_input, file = genotype_input_path, sep = "\t", row.names = FALSE, quote = FALSE)

  for (report in names(REPORTS)) {
    tmp_dir <- tempfile(paste0(report, "_retained_"))
    enchantr_report(report, report_params = c(list(
      "input" = genotype_input_path, "imgt_db" = novel_db, "species" = "human",
      "outdir" = tmp_dir, "log" = "test_allele_inference_command_log",
      "novel_allele_evidence" = evidence_path
    ), REPORTS[[report]]$params))

    genotype <- read.delim(file.path(tmp_dir, "enchantr", "genotypes", REPORTS[[report]]$output), sep = "\t")
    retained_novel <- genotype[grepl("_", genotype$genotyped_alleles, fixed = TRUE), ]
    expect_equal(nrow(retained_novel), 1, info = report)
    expect_equal(retained_novel$gene, "IGHV1-2", info = report)
    expect_equal(retained_novel$genotyped_alleles, "02_A318G,02,04", info = report)

    report_html <- paste(readLines(file.path(tmp_dir, "enchantr", "index.html"), warn = FALSE), collapse = "\n")
    expect_match(report_html, "Novel allele evidence", fixed = TRUE)
    expect_match(report_html, "IGHV1-2*02_A318G", fixed = TRUE)
    expect_false(grepl("no candidate novel alleles were retained in the current genotype call",
                       report_html, fixed = TRUE), info = report)
  }
})

test_that("the allele-based report always applies the depth-adjusted threshold", {
  skip_on_cran()
  skip_if_not_installed("piglet")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- tempfile("genotype_allele_based_")
  enchantr_report("genotype",
    report_params = list(
      "input" = input,
      "imgt_db" = IMGT_URL,
      "species" = "human",
      "method" = "allele_based",
      "outdir" = tmp_dir,
      "log" = "test_allele_based_command_log"
    )
  )

  out <- list.files(file.path(tmp_dir, "enchantr", "genotypes"), full.names = TRUE)
  genotype <- read.delim(grep("_genotype_report\\.tsv$", out, value = TRUE), sep = "\t")
  expect_equal(names(genotype)[seq_along(REPORT_COLUMNS)], REPORT_COLUMNS)
  expect_true("depth" %in% names(genotype))
  expect_length(grep("_threshold_provenance\\.tsv$", out), 1)
  expect_length(grep("_threshold_coverage\\.tsv$", out), 1)
})

test_that("single_clone_representative needs a clone id only when there are heavy chains", {
  skip_on_cran()
  heavy <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  expect_error(
    enchantr_report("genotype", report_params = list(
      "input" = heavy, "imgt_db" = IMGT_URL, "species" = "human",
      "single_clone_representative" = TRUE, "outdir" = tempfile("genotype_no_clone_"),
      "log" = "test_no_clone_command_log"
    )),
    "has heavy chains"
  )

  # a light chain only repertoire has no clones to represent: it is genotyped as is
  db <- suppressWarnings(airr::read_rearrangement(paste0(
    "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/",
    "testdata-reveal/sc5p_v2_hs_PBMC_1k_5fb_assigned.tsv")))
  light <- tempfile(fileext = ".tsv")
  suppressWarnings(airr::write_rearrangement(db[db$locus != "IGH", ], light))
  tmp_dir <- tempfile("genotype_light_only_")
  suppressWarnings(enchantr_report("genotype", report_params = list(
    "input" = light, "imgt_db" = IMGT_URL, "species" = "human",
    "single_clone_representative" = TRUE, "outdir" = tmp_dir,
    "log" = "test_light_only_command_log"
  )))
  genotype <- read.delim(list.files(file.path(tmp_dir, "enchantr", "genotypes"), full.names = TRUE)[1])
  expect_true(all(grepl("^IG[KL]", genotype$gene)))
})
