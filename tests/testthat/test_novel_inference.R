IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

#TODO: find a testing set
test_that("novel inference 1:1", {
  # Input in one file, output in 1 file.
  # All sequences have sampletype = FNA
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  input <- normalizePath(file.path("..", "data-tests", "novel_genotype", "data_to_test_novel_alleles.tsv"))
  tmp_dir <- file.path(tempdir(), "novel_inference_1_1")
  enchantr_report("novel_allele_inference",
    report_params = list(
                         "input" = input,
                         "imgt_db" = IMGT_URL,
                         "species" = "human",
                         "outdir" = tmp_dir,
                         "pos_range" = "1:318",
                         "nproc" = 1,
                         "log" = "test_allele_inference_command_log")
  )
  
  report_dir <- file.path(tmp_dir, "enchantr")
  novel_alleles <- read.delim(
    file.path(report_dir, "tables", "tigger-novel_novel_report.tsv"),
    sep = "\t"
  )
  expect_equal(nrow(novel_alleles), 6)
})

test_that("tigger_upper_bound reads the V call column it is given", {
  db <- data.frame(v_call = c("IGKV1-5*01", "IGKV1-5*01", "IGKV3-20*01", "IGKV1-5*01,IGKV1-6*01"),
                   v_germline_end = c(300, 310, 290, 312))
  expected <- tigger_upper_bound(db)
  expect_equal(unname(expected), c("IGKV3-20[*]", "IGKV1-5[*]"))
  renamed <- dplyr::rename(db, my_v = v_call)
  expect_identical(tigger_upper_bound(renamed, v_call = "my_v"), expected)
})
