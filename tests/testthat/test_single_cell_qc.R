toy_db <- data.frame(
    sample_id = c("id1", "id1", "id2", "id2", "id2", "id2"),
    sequence_id = c("seq1", "seq2", "seq3", "seq4", "seq5", "seq6"),
    cell_id   = c("cell1", "cell1", "cell2", "cell3", "cell3", "cell3"),
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
        'sample_id' = c("s1","s1","s1", "s1", "s2","s2","s2","s2", "s3", "s4"),
        'cell_id' = c("cell_1", "cell_2", "cell_3", "cell_3","cell_1", "cell_2", "cell_4","cell_5", "cell_2", "cell5"),
        'sequence_alignment' = c("AA","AA","CC", "CC", "AA","TT","CC","GG", "AM", "GGG")
    ) 
    db$sequence_id <- paste0("seq_",as.character(1:nrow(db)))
    expected_removed <- c(
        "seq_1", "seq_5",# Duplicated in samples s1 and s2
        "seq_2", "seq_9" # In s1 (AA) and s2 (AM, degenerate). Not "seq_6", because difference sequence.
        # "seq_3" and "seq_4" stay because duplicated in the same group s1
    )
    
    dups <- findSingleCellDuplicates(db, 
                             fields="sample_id", 
                             cell_id="cell_id", 
                             seq="sequence_alignment",
                             sequence_id="sequence_id")
    expect_equal(sort(dups$dups %>% filter(sc_duplicate) %>% pull(sequence_id)),
                 sort(expected_removed))
    
    db <- removeSingleCellDuplicates(db, 
                                     groups="sample_id", 
                                     cell="cell_id", 
                                     seq="sequence_alignment",
                                     sequence_id="sequence_id")
    expect_equal(length(intersect(db[['sc_duplicate']], expected_removed)), 0)
})

test_that("countSequencesPerCell", {
    countSequencesPerCell(toy_db)
})

test_that("findLightOnlyCells", {
    expect_equal(findLightOnlyCells(toy_db)[['light_only_cell']], c(F,F,T,F,F,F))
})

test_that("removeDoublets", {
    
    expected <- c("seq1", "seq2","seq3")
    expect_equal(removeDoublets(toy_db, fields='sample_id')[['sequence_id']], expected)
    
    expected <- c("seq1", "seq2", "seq1")
    expect_equal(removeDoublets(toy_db_dup_seqid, fields='sample_id')[['sequence_id']], expected)
    
    expected <- c("seq1", "seq2", "seq1")
    expect_equal(removeDoublets(toy_db_dup_cellid, fields='sample_id')[['sequence_id']], expected)
    
    # not considering sample_id
    expected <- c("seq1")
    expect_equal(removeDoublets(toy_db_dup_cellid)[['sequence_id']], expected)
})