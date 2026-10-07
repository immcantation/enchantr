test_that("a sequence with no clone is its own representative", {
  db <- data.frame(
    sequence_alignment = c("ACGT", "ACGA", "ACGT", "ACGA", "ACGT"),
    germline_alignment_d_mask = "ACGT",
    clone_id = c("1", "1", "", "", NA),
    locus = c("IGH", "IGH", "IGK", "IGK", "IGL"),
    stringsAsFactors = FALSE
  )

  out <- select_clonal_representatives(db)

  # the heavy clone keeps only its unmutated member
  expect_equal(out$clone_representative[out$locus == "IGH"], c(TRUE, FALSE))
  # every sequence without a clone_id survives, mutated or not
  expect_true(all(out$clone_representative[out$locus != "IGH"]))
  expect_equal(out$clone_size[out$locus != "IGH"], rep(1L, 3))
  expect_false("unclonal" %in% names(out))
})
