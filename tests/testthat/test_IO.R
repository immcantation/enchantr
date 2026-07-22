test_that("read input url", {
    
    skip_on_cran()
    
    url <- "https://zenodo.org/records/8179846/files/BCR_data_sample1.tsv"
    # TODO: expect no errors or warning once file in Zenodo updated
    # with rev_comp = FALSE
    # Older airr versions warn that rev_comp is not logical; newer versions
    # don't. Muffle the warning if it happens, but let any other (unexpected)
    # warning propagate and fail the test.
    withCallingHandlers(
        readInput(url),
        warning = function(w) {
            if (grepl("rev_comp is not logical for row\\(s\\):", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
    
    url <- "https://zenodo.org/records/8179846/files/BCR_data.tsv"
    expect_warning(readInput(url), "sequence_id\\(s\\) are not unique")
})
