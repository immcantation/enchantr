test_that("locus_matches takes a class, one locus or several", {
    x <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG", NA)
    expect_true(locus_matches("IGI", "IG"))
    expect_true(locus_matches("IGI", "IGI"))
    expect_equal(locus_matches(x, "IG"), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
    expect_equal(locus_matches(x, "IGH"), c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
    expect_equal(locus_matches(x, "IGK, igl"), c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
    expect_equal(locus_matches(x, "TR"), c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
    expect_equal(locus_matches(x, "TRB"), c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
    expect_error(locus_matches(x, "IGX"), "Unknown locus")
    expect_error(locus_matches(x, ""), "Unknown locus")
})

# Render the report on db and return what it wrote.
run_annotation_qc <- function(db, ...) {
    tmp <- tempfile("annotation_qc_")
    dir.create(tmp)
    input <- file.path(tmp, "S1.tsv")
    suppressWarnings(airr::write_rearrangement(db, input))
    suppressWarnings(enchantr_report("annotation_qc", report_params = list(
        input = input, outname = "S1", outdir = tmp, log = "S1_qc_command_log", ...)))
    repertoires <- file.path(tmp, "enchantr", "repertoires")
    list(
        pass = suppressWarnings(airr::read_rearrangement(file.path(repertoires, "S1_qc-pass.tsv"))),
        productive_fail = file.path(repertoires, "S1_productive-fail.tsv"),
        locus = read.delim(file.path(tmp, "enchantr", "tables", "locus_summary.tsv")),
        tables = file.path(tmp, "enchantr", "tables"),
        html = paste(readLines(file.path(tmp, "enchantr", "index.html"), warn = FALSE), collapse = "\n"),
        log = readLines(file.path(tmp, "S1_qc_command_log.txt"), warn = FALSE)
    )
}

test_that("annotation_qc report applies each filter and counts it", {
    skip_on_cran()
    ok <- strrep("A", 250)
    # rows: 1 pass, 2 unproductive, 3 no call, 4 IGK, 5 locus mismatch,
    #       6 12% N, 7 too short, 8 out of frame, 9 D and J but no V call
    db <- data.frame(
        sequence_id = paste0("s", 1:9),
        v_call = c("IGHV1-2*02", "IGHV1-2*02", NA, "IGKV1-5*01", "IGHV1-2*02", "IGHV1-2*02", "IGHV1-2*02", "IGHV1-2*02", NA),
        d_call = c(NA, NA, NA, NA, NA, NA, NA, NA, "IGHD3-10*01"),
        j_call = c("IGHJ4*02", "IGHJ4*02", NA, "IGKJ1*01", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02"),
        locus = c("IGH", "IGH", "IGH", "IGK", "IGK", "IGH", "IGH", "IGH", "IGH"),
        sequence_alignment = c(ok, ok, ok, ok, ok, paste0(strrep("N", 30), strrep("A", 220)), strrep("A", 150), ok, ok),
        productive = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        junction_length = c(45, 45, 45, 33, 45, 45, 45, 44, 45),
        stringsAsFactors = FALSE
    )
    out <- run_annotation_qc(db, locus = "IGH", restrict_to_locus = TRUE)

    expect_equal(out$pass$sequence_id, "s1")
    # s3 (no call) and s9 (D and J only) go at the V call step, not as another locus
    expect_match(out$html, "2 sequences without a V call were removed")
    expect_match(out$html, "1 sequences whose V call is from another locus were removed")
    expect_equal(out$locus$input[out$locus$locus == "IGH"], 7)
    expect_equal(out$locus$output[out$locus$locus == "IGK"], 0)
    # at the locus step: s1, s2, s8 (IGH) and s4 (IGK) are left
    locus_filter <- read.delim(file.path(out$tables, "locus_filter.tsv"))
    expect_equal(locus_filter$locus, c("IGH", "IGK"))
    expect_equal(locus_filter$sequences, c(3, 1))
    expect_equal(locus_filter$retained, c(3, 0))
    expect_equal(locus_filter$targeted, c(TRUE, FALSE))
    # non-productive and out-of-frame sequences are both kept aside
    expect_equal(suppressWarnings(airr::read_rearrangement(out$productive_fail))$sequence_id, c("s2", "s8"))
    expect_equal(out$log, c("START> AnnotationQC", "FILE> S1.tsv", "OUTPUT> S1_qc-pass.tsv",
                            "PASS> 1", "FAIL> 8"))
})

test_that("annotation_qc report keeps every row when its filters are off", {
    skip_on_cran()
    db <- data.frame(sequence_id = "s1", v_call = "IGKV1-5*01", d_call = NA, j_call = "IGKJ1*01",
                     sequence_alignment = "AC", productive = FALSE, junction_length = 1)
    out <- run_annotation_qc(db, skip_alignment_filter = TRUE, productive_only = FALSE)
    expect_equal(nrow(out$pass), 1)
    expect_equal(out$pass$locus, "IGK")
    expect_false(file.exists(out$productive_fail))
})
