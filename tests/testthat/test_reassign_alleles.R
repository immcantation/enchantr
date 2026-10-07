IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

test_that("reassign_alleles_after_novel", {
  # Input in 3 files, output in 4 file.
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

  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)
  db <- read_rearrangement(repertoires)
  expect_equal(length(grep("_", db$v_call)), 1259)

  # With the novel allele table instead of the reference, only reads carrying a
  # novel allele's parent may change
  novel_table <- list.files(file.path(novel_report_dir, "tables"), "_novel_report\\.tsv$", full.names = TRUE)
  tmp_dir <- file.path(tempdir(), "genotype_novel_reassign_table")
  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = input,
      "imgt_db" = novel_table,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "segments" = "v",
      "log" = "test_reassign_alleles_table_command_log"
    )
  )
  db_in <- read_rearrangement(input)
  db_out <- read_rearrangement(list.files(file.path(tmp_dir, "enchantr", "repertoires"), full.names = TRUE))
  expect_equal(nrow(db_out), nrow(db_in))
  before <- db_in$v_call[match(db_out$sequence_id, db_in$sequence_id)]
  changed <- before != db_out$v_call
  parents <- read.delim(novel_table)$germline_call
  expect_true(any(grepl("_", db_out$v_call)))
  expect_true(all(vapply(strsplit(before[changed], ","), function(x) any(x %in% parents), logical(1))))
})


test_that("reassign_alleles_after_genotype_inference", {
  # Input in 3 files, output in 4 file.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files", "db_let_12.tsv"))
  tmp_dir <- file.path(tempdir(), "tigger_bayesian_genotype_reassign")
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
  tmp_dir <- file.path(tempdir(), "reassign_alleles_after_genotype_inference")
  genotype_db <- file.path(report_dir, "references", "sample", "db_genotype")

  enchantr_report("reassign_alleles",
    report_params = list(
      "input" = input,
      "imgt_db" = genotype_db,
      "species" = "human",
      "outputby" = "subject_id",
      "outdir" = tmp_dir,
      "log" = "test_reassign_alleles_command_log"
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)
  db <- read_rearrangement(repertoires)
  expect_equal(unique(grep("IGHJ3", db$j_call, value = TRUE)), "IGHJ3*02")
  expect_equal(sum(db$j_call == "IGHJ3*02"), 172)
})
