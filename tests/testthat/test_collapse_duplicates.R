#### collapseDuplicates function ####

# test sequences were group by sample_id, v_call, d_call, j_call, junction_length, productive and seq_len by default
# test num_files annotation consensus_count is summed
test_that("findDuplicates", {
  # Example data.frame
  db <- data.frame(sequence_id=LETTERS[1:9],
                   sequence_alignment=c("CCCCTGGG", "CCCCTGGN","ATCGGGGA", "NAACTGGN", "NNNCTGNN", "NAACTGNG","CCCCTGGG", "CCCCTGGN","ATCGGGGA"),
                   sample_id=c('S1','S1','S1','S2','S2','S2','S2','S2','S2'),
                   v_call=c(rep("IGHV1",6),rep("IGHV2",3)),
                   d_call=rep("d_call",9),
                   j_call=rep("IGHJ2",9),
                   junction_length=rep(48,9), 
                   productive=rep(TRUE,9),
                   consensus_count=seq(1:9),
                   stringsAsFactors=FALSE)
  
  obs <- findDuplicates(db,num_fields=c('consensus_count'))
  
  exp <- data.frame(
    seq_len=rep(8,9),
    sample_id=c('S1','S1','S1','S2','S2','S2','S2','S2','S2'),
    sequence_id=LETTERS[1:9],
    sequence_alignment=c("CCCCTGGG", "CCCCTGGN","ATCGGGGA", "NAACTGGN", "NNNCTGNN", "NAACTGNG","CCCCTGGG", "CCCCTGGN","ATCGGGGA"),
    v_call=c(rep("IGHV1",6),rep("IGHV2",3)),
    d_call=rep("d_call",9),
    j_call=rep("IGHJ2",9),
    junction_length=rep(48,9),
    productive=rep(TRUE,9),
    collapse_pass=c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE),
    collapse_count=c(2,NA,1,3,NA,NA,2,NA,1),
    consensus_count=c(3,NA,3,15,NA,NA,15,NA,9),
    stringsAsFactors = F
  )
  
  expect_equivalent(obs, exp)
})


#### Collapse Duplicates Report ####
test_that("Collapse duplicates on sample_id", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample")
  report_params <- list(
    input = input,
    collapseby = "sample_id", outputby = "sample_id",
    outdir = tmp_dir,
    nproc = 1,
    log = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  print(repertoires)
  # test number and name of output files, test number of sequence and collapse_count in file S2_collapse_collapse-pass.tsv
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv','S2_collapse_collapse-pass.tsv','S3_collapse_collapse-pass.tsv','S4_collapse_collapse-pass.tsv'))
  db_s2 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S2_collapse_collapse-pass.tsv")))
  expect_equal(nrow(db_s2), 3)
  expect_equal(db_s2$collapse_count, c('3','2','1'))
})


test_that("Collapse duplicates on subject_id", {
  # Input in one file, output in 2 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir_subject <- file.path(tempdir(),"collapse_duplicates_on_subject")
  report_params <- list(
    input = input,
    collapseby = "subject_id", outputby = "subject_id",
    outdir = tmp_dir_subject,
    nproc = 1,
    log = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir_subject, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  print(repertoires)
  # test number of output files, name of output files, number of sequence , collapse_count and consensus count in file Subject_B_collapse_collapse-pass.tsv
  expect_equal(length(repertoires), 2)
  expect_equal(basename(repertoires),c('Subject_A_collapse_collapse-pass.tsv','Subject_B_collapse_collapse-pass.tsv'))
  db_B <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "Subject_B_collapse_collapse-pass.tsv")))
  expect_equal(nrow(db_B), 5)
  expect_equal(db_B$collapse_count, c('1','1','1','2','4'))
  expect_equal(db_B$consensus_count, c(10,11,12,29,64))
})


test_that("Collapse duplicates on sample_id, filter out collapse_count<2", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample_filter_count")
  report_params <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'collapse_filter_threshold' = 2,
    'outdir' = tmp_dir,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence and collapse_count in file S2_collapse_collapse-pass.tsv
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv','S2_collapse_collapse-pass.tsv','S3_collapse_collapse-pass.tsv','S4_collapse_collapse-pass.tsv'))
  db_s2 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S2_collapse_collapse-pass.tsv")))
  expect_equal(nrow(db_s2), 2)
  expect_equal(db_s2$collapse_count, c('3','2'))
})


test_that("Collapse duplicates on sample_id, mask at IMGT position 3", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir <- file.path(tempdir(),"collapse_duplicates_on_sample_IMGT_masking")
  report_params <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'mask_imgt_position' = 3,
    'outdir' = tmp_dir,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params))
  report_dir <- file.path(tmp_dir, "enchantr")
  repertoires <- list.files(file.path(report_dir, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence, collapse_count and consunsus count in S3_collapse_collapse-pass.tsv
  expect_equal(length(repertoires), 4)
  expect_equal(basename(repertoires),c('S1_collapse_collapse-pass.tsv','S2_collapse_collapse-pass.tsv','S3_collapse_collapse-pass.tsv','S4_collapse_collapse-pass.tsv'))
  db_s3 <- suppressWarnings(read_rearrangement(file.path(report_dir, "repertoires", "S3_collapse_collapse-pass.tsv")))
  expect_equal(nrow(db_s3), 2)
  expect_equal(db_s3$collapse_count, c('4','2'))
  expect_equal(db_s3$consensus_count, c(50,25))
})


test_that("Collapse duplicates on sample_id, then collapse duplicates on subject_id ", {
  # Input in one file, output in 4 files.
  skip_on_cran()
  input <- normalizePath(file.path("..", "data-tests", "subj_multiple_files","data_to_test_collapse_duplicates.tsv"))
  tmp_dir1 <- file.path(tempdir(),"collapse_duplicates_on_sample")
  report_params1 <- list(
    'input' = input,
    'collapseby' = "sample_id", 'outputby' = "sample_id",
    'outdir' = tmp_dir1,
    'nproc' = 1,
    'log' = "test_collapse_duplicate_command_log"
  )
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params1))
  
  tmp_dir2 <- file.path(tempdir(),"collapse_duplicates_on_sample_then_subject")
  report_params2 <- list(
    input = file.path(tmp_dir1, "enchantr","repertoires"),
    collapseby = "subject_id", outputby = "subject_id",
    collapse_count_colname = 'sample_count',
    outdir = tmp_dir2,
    nproc = 1,
    log = "test_collapse_duplicate_command_log"
  )  
  suppressWarnings(enchantr_report("collapse_duplicates", report_params = report_params2))
  report_dir2 <- file.path(tmp_dir2, "enchantr")
  repertoires2 <- list.files(file.path(report_dir2, "repertoires"), full.names = TRUE)
  # test number and name of output files, test number of sequence and collapse_count in file Subject_B_collapse_collapse-pass.tsv
  expect_equal(length(repertoires2), 2)
  expect_equal(basename(repertoires2),c('Subject_A_collapse_collapse-pass.tsv','Subject_B_collapse_collapse-pass.tsv'))
  db_B <- suppressWarnings(read_rearrangement(file.path(report_dir2, "repertoires", "Subject_B_collapse_collapse-pass.tsv")))
  expect_equal(nrow(db_B), 5)
  expect_equal(db_B$collapse_count, c('1','1','1','2','4'))
  expect_equal(db_B$consensus_count, c(10,11,12,29,64))
  expect_equal(db_B$sample_count, c('1','1','1','2','2'))
})