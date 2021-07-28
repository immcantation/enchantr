#' Create an Immcantation input validation project
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Input Validation
#' to create the skeleton of an Immcantation Input Validation project
#' @param  path path to the directory where the project will be created
validate_input_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "input_validation_project_files")
    project_dir <- path
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 

#' Render an Immcantation project
#' 
#' @param  name report name
#' @param  params params list, needed by report. must include outdir
enchantr_report <- function(name, report_params) {
    
    if (!dir.exists(report_params$outdir)) {
        dir.create(report_params$outdir, recursive = T)
    }
    
    outdir <- normalizePath(report_params$outdir)
    
    # Create project in outdir
    switch (name,
            "validate_input" = invisible(validate_input_project(outdir))
    )

    # render
    xfun::in_dir(
        outdir,
        book <- bookdown::render_book(
            "index.Rmd",
            output_format='bookdown::gitbook',
            config_file = "_bookdown.yml",
            clean=FALSE,
            new_session=FALSE,
            params=report_params)
    )
    book
}
    