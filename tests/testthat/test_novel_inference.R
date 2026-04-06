imgt_url <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

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
                         "imgt_db" = imgt_url,
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

  plot_paths <- list.files(
    file.path(report_dir, "ggplots"),
    pattern = "^novel-allele-[0-9]+\\.RData$",
    full.names = TRUE
  )
  expect_equal(length(plot_paths), nrow(novel_alleles))

  load(plot_paths[[1]])
  expect_true(exists("p_obj"))
  expect_true(exists("p_obj_extra_objects"))
  expect_s3_class(p_obj, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p_obj))
})
