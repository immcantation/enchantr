IMGT_URL <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"
IMGT_DB <- prepareIMGT(IMGT_URL)

# test_that("plotDbOverlap", {
#     
#     db <- data.frame(
#         'sample_id' = c("s1","s1","s1","s2","s2","s2","s2"),
#         'junction' = c("AAA","AAT","GTG","TAA","TTA","AAT","AAA"),
#         'clone_id' = c(1,1,2,3,4,5,6)
#     ) 
# 
#     
#     overlap <- plotDbOverlap(db, group="sample_id", 
#                              features=c("junction","clone_id"), 
#                              identity=c("ham_nt","exact"), 
#                              similarity=c("jaccard","jaccard"),
#                              threshold=0)
#  
# })

test_that("define clones 1:1", {
    # Input in one file, output in 1 file.
    
    skip_on_cran()
    input <- file.path("..", "data-tests", "subj_multiple_files")[1]
    input <- paste(normalizePath(list.files(input, full.names = T)),
                   collapse=",")
    tmp_dir <- file.path(tempdir(),"define_clones_1_1")
    enchantr_report('define_clones', 
                              report_params=list('input'=input, 
                                                 'imgt_db'=IMGT_DB, 
                                                 'species'='human', 
                                                 'cloneby'='sampletype', 'outputby'='sampletype',
                                                 'force'=FALSE, 
                                                 'threshold'=0.12, 
                                                 'singlecell'=NULL,
                                                 'outdir'=tmp_dir,
                                                 'nproc'=1,
                                                 'log'='test_clone_command_log'))
    report_dir <- file.path(tmp_dir,"enchantr")
    repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = T)
    expect_equal(basename(repertoires), "FNA_define-clones_clone-pass.tsv")
    db <- read_rearrangement(repertoires)
    expect_equal(nrow(db), 2998)

})

test_that("define clones 1:n", {
    # Input in one file, output in multiple files.
    
    skip_on_cran()
    input <- file.path("..", "data-tests", "subj_multiple_files")[1]
    input <- paste(normalizePath(list.files(input, full.names = T)),
                   collapse=",")
    tmp_dir <- file.path(tempdir(),"define_clones_1_n")
    enchantr_report('define_clones', 
                    report_params=list('input'=input, 
                                       'imgt_db'=IMGT_DB, 
                                       'species'='human', 
                                       'cloneby'='sampletype', 'outputby'='subject_id',
                                       'force'=FALSE, 
                                       'threshold'=0.12, 
                                       'singlecell'=NULL,
                                       'outdir'=tmp_dir,
                                       'nproc'=1,
                                       'log'='test_clone_command_log'))
})