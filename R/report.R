#' Render an Immcantation project
#' 
#' @param  name report name
#' @param  report_params params list, needed by report. must include outdir
#' @export
enchantr_report <- function(name, report_params) {
    
    if (!dir.exists(report_params$outdir)) {
        message("Creating outdir ", report_params$outdir)
        dir.create(report_params$outdir, recursive = T)
    }
    
    outdir <- normalizePath(report_params$outdir)
    
    # Create project in outdir
    switch (name,
            "validate_input" = invisible(validate_input_project(outdir)),
            "file_size" = invisible(file_size_project(outdir))
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
