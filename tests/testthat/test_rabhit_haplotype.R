test_that("rabhit_haplotype report renders", {
  skip_on_cran()
  skip_if_not_installed("rabhit")

  suppressPackageStartupMessages(library(rabhit))
  data(samples_db, HVGERM, HDGERM, HJGERM, package = "rabhit")

  tmp_dir <- tempfile("rabhit_haplotype_")
  dir.create(tmp_dir)

  db <- samples_db[samples_db$subject == "I5", ]
  db$sequence_id <- seq_len(nrow(db))
  db$productive <- TRUE

  input <- file.path(tmp_dir, "input.tsv")
  airr::write_rearrangement(db, input)

  imgt_dir <- file.path(tmp_dir, "imgt_base", "human", "vdj")
  dir.create(imgt_dir, recursive = TRUE)
  tigger::writeFasta(HVGERM, file.path(imgt_dir, "imgt_human_IGHV.fasta"))
  tigger::writeFasta(HDGERM, file.path(imgt_dir, "imgt_human_IGHD.fasta"))
  tigger::writeFasta(HJGERM, file.path(imgt_dir, "imgt_human_IGHJ.fasta"))

  enchantr_report("rabhit_haplotype",
    report_params = list(
      input = input,
      imgt_db = file.path(tmp_dir, "imgt_base"),
      species = "human",
      chain = "IGH",
      toHap_col = "v_call,d_call",
      hapBy_col = "j_call",
      hapBy = "IGHJ6",
      subject_column = "subject",
      asc_annotation = FALSE,
      run_deletions = FALSE,
      outdir = tmp_dir,
      outname = "rabhit_smoke",
      log = "test_rabhit_haplotype_command_log",
      echo = FALSE,
      cache = FALSE
    )
  )

  report_dir <- file.path(tmp_dir, "enchantr")
  expect_true(file.exists(file.path(report_dir, "index.html")))

  haplotype_file <- file.path(report_dir, "haplotypes", "rabhit_smoke_haplotype.tsv")
  expect_true(file.exists(haplotype_file))
  haplotype <- read.delim(haplotype_file, sep = "\t")
  expect_gt(nrow(haplotype), 0)
  expect_true(all(c("subject", "gene") %in% names(haplotype)))

  report_html <- paste(readLines(file.path(report_dir, "index.html"), warn = FALSE), collapse = "\n")
  expect_match(report_html, "RAbHIT haplotype inference", fixed = TRUE)
  expect_match(report_html, "Haplotype inference table", fixed = TRUE)
})
