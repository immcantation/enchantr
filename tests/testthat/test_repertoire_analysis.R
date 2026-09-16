IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"
IMGT_DB <- prepareIMGT(IMGT_URL)

test_that("repertoire analysis example", {
    # Uses the report's own bundled example data and default params
    # (cloneby="subject_id", outputby="sample_id"), which splits the
    # output into one file per sample_id.
    skip_on_cran()
    tmp_dir <- file.path(tempdir(),"clonal_assignment_example")
    enchantr_report('clonal_assignment',
                    report_params=list('outdir'=tmp_dir,
                                       'nproc'=1,
                                       'log'='test_clone_command_log'))
    report_dir <- file.path(tmp_dir,"enchantr")
    repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)

    for (fn in repertoires) {
        # Load just the first 2 rows to check that the column exists
        db <- read_rearrangement(fn, n_max=2)
        expect_true("clone_id" %in% colnames(db))
    }

    tmp_dir <- file.path(tempdir(),"repertoire_analysis_example")
    enchantr_report('repertoire_analysis',
                    report_params=list('input'=repertoires[1],
                                       'outdir'=tmp_dir,
                                       'nproc'=1,
                                       'log'='test_repertoire_command_log'))
    report_repertoire_dir <- file.path(tmp_dir,"enchantr")
    repertoire_analyzed <- list.files(file.path(report_repertoire_dir, "repertoires"), full.names = T)
    for (fn in repertoire_analyzed) {
        # Load just the first 2 rows to check that the column exists
        db <- read_rearrangement(fn, n_max=2)
        expect_true("mu_freq" %in% colnames(db))
    }
})