#' Render an Immcantation project
#' 
#' @param  name report name
#' @param  report_params params list, needed by report. must include outdir
#' @export
enchantr_report <- function(name=c("validate_input", "file_size", "chimera_analysis", "single_cell_qc"), report_params=list()) {
    
    name <- match.arg(name)
    
    if (is.null(report_params[['outdir']])) {
        report_params[['outdir']] <- getwd()
    }
    if (!dir.exists(report_params[['outdir']])) {
        message("Creating outdir ", report_params[['outdir']])
        dir.create(report_params[['outdir']], recursive = T)
    }
    
    outdir <- normalizePath(report_params[['outdir']])
    report_params[['outdir']] <- outdir
    
    # Create project in outdir
    switch (name,
            "validate_input" = invisible(validate_input_project(outdir)),
            "file_size" = invisible(file_size_project(outdir)),
            "chimera_analysis" = invisible(chimera_analysis_project(outdir)),
            "single_cell_qc" = invisible(single_cell_qc_project(outdir))
    )
    
    # render
    xfun::in_dir(
        outdir,
        book <- bookdown::render_book(
            "index.Rmd",
            output_format='enchantr::immcantation',
            config_file = "_bookdown.yml",
            clean=FALSE,
            new_session=FALSE,
            params=report_params)
    )
    book
}
