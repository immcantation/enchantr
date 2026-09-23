IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

# Every method must produce this, whatever its native output looks like.
GENOTYPE_CONTRACT <- c(
  "gene", "genotyped_alleles", "sequence_count", "clone_count",
  "confidence", "note", "candidate_alleles", "counts", "total"
)

# The presentation order augment_genotype_table() writes to the report and TSV.
# `method` and `count_basis` are added there, not by infer_genotype().
REPORT_COLUMNS <- c(
  "gene", "genotyped_alleles", "sequence_count", "clone_count",
  "confidence", "note", "candidate_alleles", "counts", "total"
)

test_that("every genotype method normalizes onto the shared schema", {
  skip_on_cran()
  db <- tigger::AIRRDb
  db <- db[!is.na(db$v_call), ]

  for (method in genotype_methods()) {
    skip_if_not_installed(enchantr:::.genotype_method(method)$package)

    genotype <- infer_genotype(db, seg = "v", loci = "IGH", method = method, find_unmutated = FALSE)

    # `fraction` has no confidence measure and so emits no such column.
    required <- setdiff(GENOTYPE_CONTRACT, c("sequence_count", "clone_count",
                                             if (method == "fraction") "confidence"))
    expect_true(all(required %in% names(genotype)), info = method)
    if (method == "allele_based") {
      # PIgLET scores alleles, not genes: one z_score per genotyped allele, in order.
      expect_equal(
        lengths(strsplit(genotype$confidence, ",")),
        lengths(strsplit(genotype$genotyped_alleles, ","))
      )
    } else if (method == "bayesian") {
      # The confidence plot needs a numeric scale; TIgGER reports k_diff as character.
      expect_type(genotype$confidence, "double")
    } else {
      expect_null(genotype$confidence)
    }
    expect_false("alleles" %in% names(genotype), info = method)
    # one row per gene
    expect_false(any(duplicated(genotype$gene)), info = method)
    expect_true(all(nzchar(genotype$genotyped_alleles)), info = method)
  }

  # methods should agree on which genes are present
  gene_sets <- lapply(genotype_methods(), function(m) {
    sort(infer_genotype(db, seg = "v", loci = "IGH", method = m, find_unmutated = FALSE)$gene)
  })
  expect_equal(gene_sets[[2]], gene_sets[[1]])
  expect_equal(gene_sets[[3]], gene_sets[[1]])
})

test_that("fraction has no confidence measure, bayesian and allele_based do", {
  skip_on_cran()
  db <- tigger::AIRRDb
  db <- db[!is.na(db$v_call), ]

  expect_false("confidence" %in% names(
    infer_genotype(db, seg = "v", loci = "IGH", method = "fraction", find_unmutated = FALSE)
  ))
  for (method in c("bayesian", "allele_based")) {
    skip_if_not_installed(enchantr:::.genotype_method(method)$package)
    conf <- infer_genotype(db, seg = "v", loci = "IGH", method = method, find_unmutated = FALSE)$confidence
    expect_false(any(is.na(conf)), info = method)
    expect_true(all(nzchar(conf)), info = method)
  }
})

test_that("unknown genotype method is rejected", {
  expect_error(infer_genotype(data.frame(), seg = "v", loci = "IGH", method = "nope"), "Unknown genotype method")
})

# TODO: find a testing set
test_that("genotype inference 1:1", {
  # Input in one file, output in 1 file.
  # All sequences have sampletype = FNA
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- tempfile("genotype_inference_1_1_")
  enchantr_report("genotype",
    report_params = list(
      "input" = input,
      "imgt_db" = IMGT_URL,
      "species" = "human",
      "method" = "bayesian",
      "outdir" = tmp_dir,
      "log" = "test_allele_inference_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  genotypes <- list.files(file.path(report_dir, "genotypes"), full.names = TRUE)
  genotype <- read.delim(genotypes, sep = "\t")
  expect_equal(nrow(genotype), 41)
  expect_false("alleles" %in% names(genotype))

  expected_cols <- c(REPORT_COLUMNS, "kh", "kd", "kt", "kq")
  expect_equal(names(genotype)[seq_along(expected_cols)], expected_cols)
  expect_true(all(is.na(genotype$clone_count)))
})

test_that("the deprecated genotype report names still render", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- tempfile("genotype_backcompat_")

  # tigger_bayesian_genotype must still resolve to the unified report, method=bayesian.
  expect_message(
    enchantr_report("tigger_bayesian_genotype",
      report_params = list(
        "input" = input,
        "imgt_db" = IMGT_URL,
        "species" = "human",
        "outdir" = tmp_dir,
        "log" = "test_allele_inference_command_log"
      )
    ),
    "deprecated"
  )

  genotype <- read.delim(
    list.files(file.path(tmp_dir, "enchantr", "genotypes"), full.names = TRUE),
    sep = "\t"
  )
  expect_equal(nrow(genotype), 41)
  expect_true("confidence" %in% names(genotype))
})

test_that("genotype report with novel alleles", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))

  tmp_dir <- file.path(tempdir(), "genotype_novel_inference")
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

  tmp_dir <- file.path(tempdir(), "genotype_novel_reassign")
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

  reassigned_input <- list.files(
    file.path(tmp_dir, "enchantr", "repertoires"),
    full.names = TRUE
  )
  expect_length(reassigned_input, 1)

  db <- read.delim(reassigned_input, sep = "\t", quote = "")
  target_call <- "IGHV1-2*02_A318G"
  target_gene <- sub("\\*.*", "", target_call)
  gene_db <- db[grepl(paste0("^", target_gene, "\\*"), db$v_call), ]
  target_db <- db[db$v_call == target_call, ]
  expect_gt(nrow(target_db), 0)

  genotype_input <- rbind(
    gene_db,
    target_db[rep(seq_len(nrow(target_db)), 20), ]
  )
  genotype_input$sequence_id <- paste0("genotype_novel_", seq_len(nrow(genotype_input)))
  genotype_input_path <- file.path(tempdir(), "genotype_input_with_novel.tsv")
  write.table(
    genotype_input,
    file = genotype_input_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  tmp_dir <- file.path(tempdir(), "genotype_retained_report")
  enchantr_report("genotype",
    report_params = list(
      "input" = genotype_input_path,
      "imgt_db" = novel_db,
      "species" = "human",
      "method" = "bayesian",
      "outdir" = tmp_dir,
      "log" = "test_allele_inference_command_log",
      "novel_allele_evidence" = evidence_path
    )
  )

  genotype <- read.delim(
    file.path(tmp_dir, "enchantr", "genotypes", "sample_genotype_report.tsv"),
    sep = "\t"
  )
  retained_novel <- genotype[grepl("_", genotype$genotyped_alleles, fixed = TRUE), ]
  expect_equal(nrow(retained_novel), 1)
  expect_equal(retained_novel$gene, "IGHV1-2")
  expect_equal(retained_novel$genotyped_alleles, "02_A318G,02,04")

  report_html <- paste(
    readLines(file.path(tmp_dir, "enchantr", "index.html"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(report_html, "Novel allele evidence", fixed = TRUE)
  expect_match(report_html, "IGHV1-2*02_A318G", fixed = TRUE)
  expect_false(grepl(
    "no candidate novel alleles were retained in the current genotype call",
    report_html,
    fixed = TRUE
  ))
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
