test_that("read input url", {
    
    skip_on_cran()
    
    url <- "https://zenodo.org/records/8179846/files/BCR_data_sample1.tsv"
    # TOD: expect no errors or warning once file in Zenodo updated
    # with rev_comp = FALSE
    expect_warning(readInput(url), "rev_comp is not logical for row\\(s\\):")
    #expect_no_error(readInput(url))
    
    url <- "https://zenodo.org/records/8179846/files/BCR_data.tsv"
    expect_warning(readInput(url), "sequence_id\\(s\\) are not unique")
})
