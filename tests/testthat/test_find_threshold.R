test_that("find threshold facets the distToNearest plot by cloneby", {
    # Uses the report's own bundled example data with a non-default cloneby.
    # The example data has several days per subject, so day groups are
    # distinct from the default sample_id grouping.
    skip_on_cran()
    tmp_dir <- file.path(tempdir(),"find_threshold_cloneby")
    enchantr_report('find_threshold',
                    report_params=list('cloneby'='day',
                                       'outdir'=tmp_dir,
                                       'nproc'=1,
                                       'log'='test_find_threshold_command_log'))
    report_dir <- file.path(tmp_dir,"enchantr")
    load(file.path(report_dir, "ggplots", "dtnplot.RData"))
    expect_setequal(names(dtnplot$facet$params$facets),
                    c("single_cell", "sample_id", "day"))
})
