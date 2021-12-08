test_that("file_size report runs", {
    # W: AT
    # M: AC
    db <- data.frame(
        'sample_id' = c("s1","s1","s1", "s1", "s2","s2","s2","s2", "s3"),
        'cell_id' = c("cell_1", "cell_2", "cell_3", "cell_3","cell_1", "cell_2", "cell_4","cell_5", "cell_2"),
        'sequence_alignment' = c("AA","AA","CC", "CC", "AA","TT","CC","GG", "AM")
    ) 
    db$sequence_id <- paste0("seq_",as.character(1:nrow(db)))
    expected_removed <- c(
        "seq_1", "seq_5",# Duplicated in samples s1 and s2
        "seq_2", "seq_9" # In s1 (AA) and s2 (AM, degenerate). Not "seq_6", because difference sequence.
        # "seq_3" and "seq_4" stay because duplicated in the same group s1
    )
    
    db <- findSingleCellDuplicates(db, 
                             groups="sample_id", 
                             cell="cell_id", 
                             seq="sequence_alignment",
                             sequence_id="sequence_id")
    expect_equal(sum(db[['sc_duplicate']]),length(expected_removed))
    
    db <- removeSingleCellDuplicates(db, 
                                     groups="sample_id", 
                                     cell="cell_id", 
                                     seq="sequence_alignment",
                                     sequence_id="sequence_id")
    expect_equal(length(intersect(db[['sc_duplicate']], expected_removed)), 0)
})

