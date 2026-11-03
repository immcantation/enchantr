toy_db <- data.frame(
    sample_id = c("id1", "id1", "id2", "id2", "id2", "id2"),
    sequence_id = c("seq1", "seq2", "seq3", "seq4", "seq5", "seq6"),
    cell_id = c("cell1", "cell1", "cell2", "cell3", "cell3", "cell3"),
    locus = c("IGH", "IGK", "IGK", "IGH", "IGH", "IGL")
)

toy_db_dup_seqid <- data.frame(
    sample_id = c("id1", "id1", "id2", "id2", "id2", "id2"),
    sequence_id = c("seq1", "seq2", "seq1", "seq2", "seq3", "seq4"),
    cell_id   = c("cell1", "cell1", "cell2", "cell3", "cell3", "cell3"),
    locus = c("IGH", "IGK", "IGK", "IGH", "IGH", "IGL")
)

toy_db_dup_cellid <- data.frame(
    sample_id = c("id1", "id1", "id2", "id2", "id2", "id2"),
    sequence_id = c("seq1", "seq2", "seq1", "seq2", "seq3", "seq4"),
    cell_id   = c("cell1", "cell1", "cell2", "cell1", "cell1", "cell1"),
    locus = c("IGH", "IGK", "IGK", "IGH", "IGH", "IGL")
)

test_that("sc duplicates", {
    # W: AT
    # M: AC
    db <- data.frame(
        "sample_id" = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2", "s3",
                        "s3", "s4"),
        "cell_id" = c("cell_1", "cell_2", "cell_3", "cell_3", "cell_1",
                      "cell_2", "cell_4", "cell_5", "cell_2", "cell_2",
                      "cell_5"),
        "sequence_alignment" = c("AA", "AA", "CC", "CC", "AA", "TT", "CC", "GG",
                                 "AM", "CC", "GGG")
        )
    db$sequence_id <- paste0("seq_",as.character(1:nrow(db)))
    db$locus <- "IGH"
    dup_seqs <- c(
        "seq_1", "seq_5",# Duplicated in samples s1 and s2
        "seq_2", "seq_9" # In s1 (AA) and s2 (AM, degenerate). Not "seq_6", because difference sequence.
        # "seq_3" and "seq_4" stay because duplicated in the same group s1
    )

    # seq10 is in a cell with another sequence that is found in other samples,
    # therefore seq10 will be removed using the 'cells' mode in removeSingleCellDuplicates
    dup_cells <- c("seq1", "seq2", "seq5", "seq9", "seq10")

    dups <- findSingleCellDuplicates(db,
                                     fields = "sample_id",
                                     cell_id = "cell_id",
                                     seq = "sequence_alignment",
                                     sequence_id = "sequence_id",
                                     mode = "sequences")
    expect_equal(sort(dups$dups %>%
                          filter(sc_duplicate) %>%
                          dplyr::pull(sequence_id)),
                 sort(dup_seqs))


    db_dedup_seqs <- removeSingleCellDuplicates(db,
                                                fields = "sample_id",
                                                cell_id = "cell_id",
                                                seq = "sequence_alignment",
                                                sequence_id = "sequence_id",
                                                mode = "sequences")
    expect_equal(length(intersect(db_dedup_seqs[["sequence_id"]], dup_seqs)), 0)

    db_dedup_cells <- removeSingleCellDuplicates(db,
                                                 fields = "sample_id",
                                                 cell = "cell_id",
                                                 seq = "sequence_alignment",
                                                 sequence_id = "sequence_id",
                                                 mode = "cells")
    expect_equal(length(intersect(db_dedup_cells[["cell_id"]], dup_cells)), 0)
})

# test_that("countSequencesPerCell", {
#     countSequencesPerCell(toy_db)
# })

test_that("findLightOnlyCells", {
    expect_equal(findLightOnlyCells(toy_db)[["light_only_cell"]], c(F,F,T,F,F,F))
})

test_that("removeDoublets", {
    expected <- c("seq1", "seq2","seq3")
    expect_equal(removeDoublets(toy_db,
                                fields = "sample_id")[["sequence_id"]],
                 expected)

    expected <- c("seq1", "seq2", "seq1")
    expect_equal(removeDoublets(toy_db_dup_seqid,
                                fields = "sample_id")[["sequence_id"]],
                 expected)

    expected <- c("seq1", "seq2", "seq1")
    expect_equal(removeDoublets(toy_db_dup_cellid,
                                fields = "sample_id")[["sequence_id"]],
                 expected)

    # not considering sample_id
    expected <- c("seq1")
    expect_equal(removeDoublets(toy_db_dup_cellid)[["sequence_id"]], expected)
})


test_that("single cell qc example", {
  # Uses the report's own bundled example data and default params
  skip_on_cran()
  tmp_dir <- file.path(tempdir(),"scqc_report_example")
  enchantr_report('single_cell_qc',
                  report_params=list('outdir'=tmp_dir,
                                     'log'='test_scqc_command_log'))
  report_scqc_dir <- file.path(tmp_dir,"enchantr")
  repertoire_analyzed <- list.files(file.path(report_scqc_dir, "repertoires"), full.names = TRUE)
  expect_equal(length(repertoire_analyzed),2)
})

test_that("single cell qc contamination checking and removal are independent", {
  skip_on_cran()

  example_file <- system.file(
    "rstudio", "templates", "project", "single_cell_qc_project_files",
    "scExampleDb.tsv.gz",
    package = "enchantr"
  )
  example_db <- airr::read_rearrangement(example_file)
  example_db[["sample_id"]] <- "sample_1"

  # Copy several complete cells into a second sample. Their cell IDs and
  # sequences deliberately match sample_1, so the contamination check must
  # classify those cells as duplicates while leaving other sample_1 cells.
  eligible_db <- findLightOnlyCells(
    example_db,
    sample_id = "sample_id",
    cell_id = "cell_id",
    locus = "locus",
    fields = NULL
  )
  eligible_db <- eligible_db[!eligible_db[["light_only_cell"]], , drop = FALSE]
  eligible_db <- removeDoublets(
    eligible_db,
    cell_id = "cell_id",
    locus = "locus",
    sequence_id = "sequence_id",
    fields = "sample_id"
  )
  contaminated_cell_ids <- head(unique(eligible_db[["cell_id"]]), 10)
  expect_gt(length(contaminated_cell_ids), 0)
  contaminated_db <- example_db[
    example_db[["cell_id"]] %in% contaminated_cell_ids,
    ,
    drop = FALSE
  ]
  contaminated_db[["sample_id"]] <- "sample_2"
  test_db <- rbind(example_db, contaminated_db)
  test_input <- tempfile("scqc_contamination_input_", fileext = ".tsv.gz")
  airr::write_rearrangement(test_db, test_input)

  render_scqc <- function(checkcontamination, removecontamination) {
    tmp_dir <- tempfile("scqc_contamination_")
    dir.create(tmp_dir)

    enchantr_report(
      "single_cell_qc",
      report_params = list(
        outdir = tmp_dir,
        input = test_input,
        log = "test_scqc_command_log",
        checkcontamination = checkcontamination,
        removecontamination = removecontamination
      )
    )

    repertoire_dir <- file.path(tmp_dir, "enchantr", "repertoires")
    pass_files <- list.files(
      repertoire_dir,
      pattern = "_scqc-pass\\.tsv\\.gz$",
      full.names = TRUE
    )
    fail_files <- list.files(
      repertoire_dir,
      pattern = "_scqc-contamination-fail\\.tsv\\.gz$",
      full.names = TRUE
    )

    read_sequence_ids <- function(files) {
      if (length(files) == 0) {
        return(character())
      }
      unlist(lapply(files, function(file) {
        airr::read_rearrangement(file)[["sequence_id"]]
      }), use.names = FALSE)
    }

    list(
      pass = read_sequence_ids(pass_files),
      fail = read_sequence_ids(fail_files)
    )
  }

  unchecked <- render_scqc(
    checkcontamination = FALSE,
    removecontamination = TRUE
  )
  checked <- render_scqc(
    checkcontamination = TRUE,
    removecontamination = FALSE
  )
  removed <- render_scqc(
    checkcontamination = TRUE,
    removecontamination = TRUE
  )

  # Skipping the check retains every sequence and creates no failure output.
  expect_length(unchecked$fail, 0)

  # Check-only mode retains every sequence but saves detected duplicates.
  expect_setequal(checked$pass, unchecked$pass)
  expect_gt(length(checked$fail), 0)

  # Removal mode saves the same failures and excludes them from final output.
  expect_setequal(removed$fail, checked$fail)
  expect_setequal(
    removed$pass,
    setdiff(checked$pass, checked$fail)
  )
})
