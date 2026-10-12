#' Create an Immcantation annotation QC project
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select Immcantation annotation QC
#' to create the skeleton of an Immcantation annotation QC project
#' @param  path path to the directory where the project will be created
annotation_qc_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                              "templates", "project",
                              "annotation_qc_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }

    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}


#' Match loci against the targeted ones
#'
#' @param    x      loci, e.g. the \code{locus} column.
#' @param    locus  receptor class (\code{IG}, \code{TR}) or comma-separated loci
#'                  (\code{IGH}, \code{IGK,IGL}, \code{TRB}).
#' @return   A logical vector, \code{FALSE} for \code{NA}.
#' @examples
#' locus_matches(c("IGH", "IGK", "IGL", "TRB", NA), "IGK,IGL")
#' @export
locus_matches <- function(x, locus) {
    wanted <- toupper(trimws(strsplit(locus, ",")[[1]]))
    valid <- c("IG", "IGH", "IGK", "IGL", "TR", "TRA", "TRB", "TRD", "TRG")
    if (!all(wanted %in% valid)) {
        stop("Unknown locus '", locus, "'. Use one or more of: ", paste(valid, collapse = ", "))
    }
    grepl(paste0("^(", paste(wanted, collapse = "|"), ")"), toupper(x))
}
