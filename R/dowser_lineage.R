#' Create an Immcantation ...
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
dowser_lineage_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                              "templates", "project",
                              "dowser_lineage_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}



#' segment_call <- c("IGHA*06", "IGHD*01", "IGHG2B*02", "IGHG2C*01", "IGHM*01", "IGKC*02", "IGLC1*01", "IGLC2*01", "IGLC3*01")
#' getIsotype(segment_call)
# getIsotype <- function(segment_call, first=TRUE, collapse=TRUE,
#                        strip_d=TRUE, omit_nl=FALSE, sep=",") {
#     isotype_regex <- '((IG[HLK]|TR[ABDG]).)'
#     r <- getSegment(segment_call, isotype_regex, first=first, collapse=collapse,
#                     strip_d=strip_d, omit_nl=omit_nl, sep=sep)
#     r
#     return(r)
# }


#' @export
prepareIMGT <- function(imgt_db) {
    is_url <-  grepl("^https+://", imgt_db)
    is_zip <-  grepl(".zip", imgt_db)
    bname <- strsplit(basename(imgt_db),"\\.")[[1]][1]
    tmp_dir <- tempdir()
    if ( is_url ) {
        tmp_file <- tempfile(tmpdir=tmp_dir, fileext = ".zip")
        download.file(imgt_db, destfile=tmp_file)
        imgt_db <- tmp_file
    }
    if ( is_zip ) {
        unzip(imgt_db, exdir=tmp_dir)
        imgt_db <- file.path(tmp_dir, bname)
    }
    imgt_db
}
