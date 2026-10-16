#' Create a RAbHIT haplotype project template
#'
#' Create the directory skeleton for a RAbHIT haplotype inference project from
#' the package's bundled RStudio project template.
#'
#' @param path Path to the directory where the project will be created.
#' @param ... Additional arguments (currently unused).
#' @export
rabhit_haplotype_project <- function(path, ...) {
  skeleton_dir <- file.path(
    system.file(package = "enchantr"), "rstudio",
    "templates", "project",
    "rabhit_haplotype_project_files"
  )
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}
