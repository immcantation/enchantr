imgt_url <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

#TODO: find a testing set
test_that("novel inference 1:1", {
  # Input in one file, output in 1 file.
  # All sequences have sampletype = FNA
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- file.path(tempdir(), "novel_inference_1_1")
  enchantr_report("novel_allele_inference",
    report_params = list(
                         "input" = input,
                         "imgt_db" = imgt_url,
                         "species" = "human",
                         "outdir" = tmp_dir,
                         "pos_range" = "1:318",
                         "nproc" = 1,
                         "log" = "test_allele_inference_command_log")
  )
  
  novel_alleles <- read.delim(file.path(tmp_dir, "enchantr", "tigger_novel_report.tsv"), sep = "\t")
  expect_equal(nrow(novel_alleles), 6)
})
