#' Create an Immcantation crosscontamination detection project
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Contamination Validation
#' to create the skeleton of an Immcantation contamination project
#' @param  path path to the directory where the project will be created
contamination_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "contamination_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 
