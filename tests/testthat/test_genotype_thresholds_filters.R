test_that("threshold coverage is matched on the column piglet would use", {
  thresholds <- data.frame(
    allele = c("IGHV1-2*02", "IGHV1-3*01"),
    asc_allele = c("IGHVF1-G1*01", "IGHVF1-G1*02"),
    threshold = c(1e-04, 5e-04), stringsAsFactors = FALSE
  )
  db <- data.frame(v_call = c("IGHV1-2*02", "IGHV1-2*02", "IGHV1-9*01"), stringsAsFactors = FALSE)

  cov <- allele_threshold_coverage(db, "v_call", thresholds, asc_annotation = FALSE)
  expect_equal(nrow(cov), 2)
  expect_true(cov$matched[cov$allele == "IGHV1-2*02"])
  expect_false(cov$matched[cov$allele == "IGHV1-9*01"])

  s <- summarize_threshold_coverage(cov)
  expect_equal(s$alleles_matched, 1)
  expect_equal(s$alleles_defaulted, 1)
  expect_equal(s$sequences_matched, 2)

  # IMGT-named calls against the ASC column match nothing
  cov_asc <- allele_threshold_coverage(db, "v_call", thresholds, asc_annotation = TRUE)
  expect_equal(summarize_threshold_coverage(cov_asc)$alleles_matched, 0)
})

test_that("a threshold table without the needed name column fails loudly", {
  thresholds <- data.frame(allele = "IGHV1-2*02", threshold = 1e-04, stringsAsFactors = FALSE)
  db <- data.frame(v_call = "IGHV1-2*02", stringsAsFactors = FALSE)
  expect_error(
    allele_threshold_coverage(db, "v_call", thresholds, asc_annotation = TRUE),
    "no 'asc_allele' column"
  )
})

test_that("ambiguous calls follow the inference's assignment rule", {
  thresholds <- data.frame(allele = c("IGHV1-2*02", "IGHV1-2*04"),
                           threshold = c(1e-04, 1e-04), stringsAsFactors = FALSE)
  db <- data.frame(v_call = c("IGHV1-2*02", "IGHV1-2*02,IGHV1-2*04"), stringsAsFactors = FALSE)
  expect_equal(nrow(allele_threshold_coverage(db, "v_call", thresholds, single_assignments = TRUE)), 1)
  expect_equal(nrow(allele_threshold_coverage(db, "v_call", thresholds, single_assignments = FALSE)), 2)
})

test_that("provenance records the table actually used", {
  thresholds <- data.frame(allele = "IGHV1-2*02", threshold = 1e-04, stringsAsFactors = FALSE)
  prov <- allele_threshold_provenance(thresholds, path = NULL, release = NULL)
  expect_equal(prov$source, "piglet bundled allele_threshold_table")
  expect_equal(prov$release, "unspecified")
  expect_equal(prov$n_alleles, 1L)
  expect_match(prov$md5, "^[0-9a-f]{32}$")

  f <- tempfile(fileext = ".tsv")
  write.table(thresholds, f, sep = "\t", row.names = FALSE, quote = FALSE)
  prov2 <- allele_threshold_provenance(thresholds, path = f, release = "musa-2026-09")
  expect_equal(prov2$source, f)
  expect_equal(prov2$release, "musa-2026-09")
  expect_equal(prov2$md5, prov$md5)   # same table, same checksum, whatever its source
})

test_that("row counting reports in and out", {
  steps <- count_rows_step(NULL, "productive == TRUE", 100, 90)
  steps <- count_rows_step(steps, "junction_length %% 3 == 0", 90, 88)
  expect_equal(steps$dropped, c(10, 2))
  expect_equal(steps$rows_in[2], steps$rows_out[1])
})

test_that("a threshold table is read whether it is comma or tab separated", {
  rows <- data.frame(
    allele = c("IGHV1-2*01", "IGHV1-2*02"),
    threshold = c(1e-04, 2e-03),
    stringsAsFactors = FALSE
  )
  csv <- tempfile(fileext = ".csv")
  tsv <- tempfile(fileext = ".tsv")
  utils::write.csv(rows, csv, row.names = FALSE, quote = FALSE)
  utils::write.table(rows, tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  from_csv <- read_allele_thresholds(csv)
  from_tsv <- read_allele_thresholds(tsv)
  expect_equal(from_csv, rows)
  expect_equal(from_tsv, rows)
  expect_equal(.threshold_key_column(from_csv, asc_annotation = FALSE), "allele")
})
