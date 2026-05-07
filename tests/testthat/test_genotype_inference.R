IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

# TODO: find a testing set
test_that("genotype inference 1:1", {
  # Input in one file, output in 1 file.
  # All sequences have sampletype = FNA
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- tempfile("genotype_inference_1_1_")
  enchantr_report("tigger_bayesian_genotype",
    report_params = list(
      "input" = input,
      "imgt_db" = IMGT_URL,
      "species" = "human",
      "outdir" = tmp_dir,
      "log" = "test_allele_inference_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  genotypes <- list.files(file.path(report_dir, "genotypes"), full.names = TRUE)
  genotype <- read.delim(genotypes, sep = "\t")
  expect_equal(nrow(genotype), 41)
  expect_false("alleles" %in% names(genotype))

  expected_cols <- c(
    "gene", "genotyped_alleles", "sequence_count", "clone_count",
    "k_diff", "note", "candidate_alleles", "counts",
    "total", "kh", "kd", "kt", "kq"
  )
  expect_equal(names(genotype)[seq_along(expected_cols)], expected_cols)
  expect_true(all(is.na(genotype$clone_count)))
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
  enchantr_report("tigger_bayesian_genotype",
    report_params = list(
      "input" = genotype_input_path,
      "imgt_db" = novel_db,
      "species" = "human",
      "outdir" = tmp_dir,
      "log" = "test_allele_inference_command_log",
      "novel_allele_evidence" = evidence_path
    )
  )

  genotype <- read.delim(
    file.path(tmp_dir, "enchantr", "genotypes", "tigger_genotype_report.tsv"),
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
