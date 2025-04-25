# test_that("input_validation report runs", {
#     tmp_dir <- file.path(tempdir(),"test_input_validation","")
#     enchantr_report('validate_input', report_params=list('outdir'=tmp_dir))
#     expect_true(file.exists(file.path(tmp_dir,"enchantr","index.html")))
# })


# test_that("clonal_assignment report runs", {
#     tmp_dir <- file.path(tempdir(),"test_clonal_assignment","")
#     enchantr_report('clonal_assignment', report_params=list('outdir'=tmp_dir))
#     expect_true(file.exists(file.path(tmp_dir,"enchantr","index.html")))
# })
