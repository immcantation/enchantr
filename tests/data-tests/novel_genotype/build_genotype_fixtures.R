imgt_url <- "https://raw.githubusercontent.com/nf-core/test-datasets/airrflow/database-cache/imgtdb_base.zip"

suppressPackageStartupMessages(library(enchantr))

fixture_root <- normalizePath("tests/data-tests/novel_genotype")
input <- file.path(fixture_root, "data_to_test_novel_alleles.tsv")

copy_tree <- function(src, dst) {
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  file.copy(
    list.files(src, full.names = TRUE, all.files = TRUE, no.. = TRUE),
    dst,
    recursive = TRUE,
    overwrite = TRUE
  )
}

work_novel_dir <- tempfile("genotype_fixture_novel_")
enchantr_report("novel_allele_inference",
  report_params = list(
    input = input,
    imgt_db = imgt_url,
    species = "human",
    outdir = work_novel_dir,
    pos_range = "1:318",
    nproc = 1,
    log = "test_allele_inference_command_log"
  )
)

novel_db <- file.path(work_novel_dir, "enchantr", "db_novel")
novel_evidence <- file.path(
  work_novel_dir, "enchantr", "tigger-novel_novel_allele_evidence.rda"
)

work_reassigned_dir <- tempfile("genotype_fixture_reassigned_")
enchantr_report("reassign_alleles",
  report_params = list(
    input = input,
    imgt_db = novel_db,
    species = "human",
    outputby = "subject_id",
    outdir = work_reassigned_dir,
    segments = "v",
    log = "test_reassign_alleles_command_log"
  )
)

reassigned_input <- list.files(
  file.path(work_reassigned_dir, "enchantr", "repertoires"),
  full.names = TRUE
)
stopifnot(length(reassigned_input) == 1)
reassigned_db <- read.delim(reassigned_input[[1]], sep = "\t", quote = "")

write.table(
  reassigned_db,
  file = gzfile(file.path(fixture_root, "genotype_input_fallback.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
file.copy(
  novel_evidence,
  file.path(fixture_root, "tigger-novel_novel_allele_evidence.rda"),
  overwrite = TRUE
)
copy_tree(novel_db, file.path(fixture_root, "db_novel"))

single_assignment_novel <- reassigned_db$v_call[
  !is.na(reassigned_db$v_call) &
    !grepl(",", reassigned_db$v_call) &
    grepl("_", reassigned_db$v_call, fixed = TRUE)
]
stopifnot(length(single_assignment_novel) > 0)

target_call <- names(sort(table(single_assignment_novel), decreasing = TRUE))[1]
target_gene <- sub("\\*.*", "", target_call)

gene_rows <- reassigned_db[
  grepl(paste0("^", target_gene, "\\*"), reassigned_db$v_call),
  ,
  drop = FALSE
]
target_rows <- reassigned_db[
  reassigned_db$v_call == target_call,
  ,
  drop = FALSE
]
stopifnot(nrow(target_rows) > 0)

synthetic_db <- rbind(
  gene_rows,
  target_rows[rep(seq_len(nrow(target_rows)), 20), , drop = FALSE]
)
synthetic_db$sequence_id <- paste0("synthetic_retained_novel_", seq_len(nrow(synthetic_db)))

write.table(
  synthetic_db,
  file = gzfile(file.path(fixture_root, "genotype_input_retained.tsv.gz")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
