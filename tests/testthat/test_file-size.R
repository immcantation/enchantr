test_that("file size report issue 49", {
    # https://github.com/immcantation/enchantr/issues/49
    # nextflow pull Vivian0105/airrflow -r new
    # nextflow run Vivian0105/airrflow -r new -profile test,docker \
    # --outdir test-airrflow-results \
    # -resume -w work-test-airrflow
    # tests/testthat/issue-49/work-test-airrflow/d0/2e4e38a52dc9d51f37226671055eac
    # Error was:
    #
    # 12/26 [unnamed-chunk-5] 
    # Error in `mutate()`:
    # i In argument: `name = ifelse(n > 1, .makeVertexName(.data), input)`.
    # i In row 7.
    # Caused by error in `.makeVertexName()`:
    # ! Unexpected duplicated input names.
    # Backtrace:
    #     x
    #     1. +-enchantr::consoleLogsAsGraphs(log_data, metadata)
    #     2. | \-logs %>% rowwise() %>% ...
    #     3. +-dplyr::mutate(...)
    #     4. +-dplyr:::mutate.data.frame(., name = ifelse(n > 1, .makeVertexName(.data), input))
    #     5. | \-dplyr:::mutate_cols(.data, dplyr_quosures(...), by)
    #     6. |   +-base::withCallingHandlers(...)
    #     7. |   \-dplyr:::mutate_col(dots[[i]], data, mask, new_columns)
    #     8. |     \-mask$eval_all_mutate(quo)
    #     9. |       \-dplyr (local) eval()
    # 10. +-base::ifelse(n > 1, .makeVertexName(.data), input)
    # 11. +-enchantr (local) .makeVertexName(.data)
    # 12. | \-base::stop("Unexpected duplicated input names.")
    # 13. \-base::.handleSimpleError(...)
    # 14.   \-dplyr (local) h(simpleError(msg, call))
    # 15.     \-rlang::abort(message, class = error_class, parent = parent, call = error_call)
    
    log_table <- structure(
        list(step = c(1, 1, 1, 1, 1, 1, 1, 1, 1), 
             task = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
                              levels = "ParseDb-split", class = "factor"), 
             field = c("FILE", "FIELD", "NUM_SPLIT", "OUTPUT1", "OUTPUT2", 
                       "RECORDS", "PARTS", "SIZE1", "SIZE2"), 
             value = c("Sample8_quality-pass.tsv", "productive", "None", 
                       "Sample8_productive-F.tsv", "Sample8_productive-T.tsv", 
                       "92", "2", "3", "89")), 
        class = "data.frame", 
        row.names = 3:11)
    out_log <-     formatConsoleLog(log_table)
    expeted_log <- structure(
        list(
            input = c("Sample8_quality-pass.tsv", "Sample8_quality-pass.tsv"), 
            output = c("Sample8_productive-F.tsv", "Sample8_productive-T.tsv"),
            task = structure(c(1L, 1L), class = "factor", levels = "ParseDb-split"), 
            input_size = c(92, 92), 
            output_size = c(3, 89)), 
            row.names = c(NA, -2L), 
            class = "data.frame")
    expect_equal(out_log, expeted_log)

})

test_that("makedb igblast parses OUTPUT>", {

    log_text <- c("         START> MakeDB", "       COMMAND> igblast", "  ALIGNER_FILE> BLOD_AM1_2_A1_igblast.fmt7", 
                    "      SEQ_FILE> BLOD_AM1_2_A1_parse-select_sequences.fasta", 
                    "       ASIS_ID> False", "    ASIS_CALLS> False", "      VALIDATE> strict", 
                    "      EXTENDED> True", "INFER_JUNCTION> False", "PROGRESS> 09:59:31 |Loading files       | 0.0 min", 
                    "PROGRESS> 09:59:35 |Done                | 0.1 min", "PROGRESS> 09:59:35 |                    |   0% (      0) 0.0 min", 
                    "PROGRESS> 09:59:43 |#                   |   5% (  9,013) 0.1 minOUTPUT> BLOD_AM1_2_A1_db-pass.tsv", 
                    "  PASS> 13353", "  FAIL> 38", "   END> MakeDb")
    out_log <- loadConsoleLog(log_text)
    expected_log <- data.frame(
        step = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        task = c("MakeDB-igblast", "MakeDB-igblast", "MakeDB-igblast", "MakeDB-igblast", 
                 "MakeDB-igblast", "MakeDB-igblast", "MakeDB-igblast", "MakeDB-igblast", 
                 "MakeDB-igblast", "MakeDB-igblast"),
        field = c("ALIGNER_FILE", "SEQ_FILE", "ASIS_ID", "ASIS_CALLS", "VALIDATE", 
                  "EXTENDED", "INFER_JUNCTION", "PASS", "FAIL", "OUTPUT"),
        value = c("BLOD_AM1_2_A1_igblast.fmt7", "BLOD_AM1_2_A1_parse-select_sequences.fasta",
                  "False", "False", "strict", "True", "False", "13353", "38", 
                  "BLOD_AM1_2_A1_db-pass.tsv"),
        stringsAsFactors = FALSE
    )
    expect_equivalent(as.matrix(out_log), as.matrix(expected_log))
    
})
