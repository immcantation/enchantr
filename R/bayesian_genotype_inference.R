# The tigger_bayesian_genotype report, kept for callers that use it (e.g.
# nf-core/airrflow). It runs infer_genotype(method = "bayesian").

#' Create a TIGGER Bayesian genotype project template
#'
#' Create the directory skeleton for a TIGGER Bayesian genotype project from
#' the package's bundled RStudio project template. This is a convenience
#' wrapper used by RStudio's New Project > New Directory > From Package
#' templates.
#'
#' @param path Path to the directory where the project will be created.
#' @param ... Additional arguments (currently unused).
#' @export
tigger_bayesian_genotype_project <- function(path, ...) {
  skeleton_dir <- file.path(
    system.file(package = "enchantr"), "rstudio",
    "templates", "project",
    "tigger_bayesian_genotype_project_files"
  )
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}
