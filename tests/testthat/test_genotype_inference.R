imgt_url <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

#TODO: find a testing set
test_that("novel inference 1:1", {
  # Input in one file, output in 1 file.
  # All sequences have sampletype = FNA
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- file.path(tempdir(), "genotype_inference_1_1")
  enchantr_report("tigger_bayesian_genotype",
    report_params = list(
                         "input" = input,
                         "imgt_db" = imgt_url,
                         "species" = "human",
                         "outdir" = tmp_dir,
                         "call_column" = "v_call",
                         "reassign" = TRUE,
                         "log" = "test_allele_inference_command_log")
  )
  
  genotype <- read.delim(file.path(tmp_dir, "enchantr", "tigger_v_call_genotype_report.tsv"), sep = "\t")
  expect_equal(nrow(genotype), 10)
})
