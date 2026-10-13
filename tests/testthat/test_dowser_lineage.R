test_that("dowser lineage example keeps all default num_fields", {
    # Uses the report's own bundled example data and default params.
    # build="pratchet" so that the test doesn't need the IgPhyML executable.
    skip_on_cran()
    tmp_dir <- file.path(tempdir(),"dowser_lineage_example")
    suppressWarnings(
        enchantr_report('dowser_lineage',
                        report_params=list('outdir'=tmp_dir,
                                           'build'='pratchet',
                                           'nproc'=1,
                                           'log'='test_dowser_command_log'))
    )
    report_dir <- file.path(tmp_dir,"enchantr")

    clones <- readRDS(file.path(report_dir, "formatted_clones.RDS"))
    expect_gt(nrow(clones), 0)
    # Both default num_fields must reach the clone data. consensus_count was
    # silently dropped when the default had a space after the comma.
    for (i in seq_len(nrow(clones))) {
        clone_fields <- colnames(clones[['data']][[i]]@data)
        expect_true(all(c("duplicate_count", "consensus_count") %in% clone_fields))
    }

    trees <- readRDS(file.path(report_dir, "trees_object.RDS"))
    expect_equal(nrow(trees), nrow(clones))
})
