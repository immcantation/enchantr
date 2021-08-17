test_that("input_validation report runs", {
    tmp_dir <- file.path(tempdir(),"test_input_validation","")
    enchantr_report('validate_input', report_params=list('outdir'=tmp_dir))
})
