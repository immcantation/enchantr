#' Render an Immcantation project
#'
#' @param  name report name
#' @param  report_params params list, needed by report. must include outdir
#' @param  help if \code{TRUE}, print the parameters accepted by report
#'               \code{name}, with their default values and descriptions,
#'               and return them invisibly without rendering the report.
#'
#' @examples
#' enchantr_report("collapse_duplicates", help=TRUE)
#' @export
enchantr_report <- function(name=c("validate_input",
                                   "file_size",
                                   "chimera_analysis",
                                   "single_cell_qc",
                                   "annotation_qc",
                                   "contamination",
                                   "collapse_duplicates",
                                   "find_threshold",
                                   "clonal_assignment",
                                   "repertoire_analysis",
                                   "convergence",
                                   "dowser_lineage",
                                   "novel_allele_inference",
                                   "genotype",
                                   "rabhit_haplotype",
                                   "tigger_bayesian_genotype",
                                   "reassign_alleles"), report_params=list(),
                            help=FALSE) {

    name <- match.arg(name)

    if (help) {
        return(invisible(print_report_params(name)))
    }

    if (is.null(report_params[['outdir']])) {
        report_params[['outdir']] <- getwd()
    }
    if (!dir.exists(report_params[['outdir']])) {
        message("Creating outdir ", report_params[['outdir']])
        dir.create(report_params[['outdir']], recursive = T)
    }
    
    # outdir <- normalizePath(report_params[['outdir']])
    outdir <- report_params[['outdir']]
    ##report_params[['outdir']] <- outdir
    
    # Create project in outdir
    switch (name,
            "validate_input" = invisible(validate_input_project(outdir)),
            "file_size" = invisible(file_size_project(outdir)),
            "chimera_analysis" = invisible(chimera_analysis_project(outdir)),
            "single_cell_qc" = invisible(single_cell_qc_project(outdir)),
            "annotation_qc" = invisible(annotation_qc_project(outdir)),
            "contamination" = invisible(contamination_project(outdir)),
            "collapse_duplicates" = invisible(collapse_duplicates_project(outdir)),
            "find_threshold" = invisible(find_threshold_project(outdir)),
            "clonal_assignment" = invisible(clonal_assignment_project(outdir)),
            "repertoire_analysis" = invisible(repertoire_analysis_project(outdir)),
            "convergence" = invisible(convergence_project(outdir)),            
            "dowser_lineage" = invisible(dowser_lineage_project(outdir)),
            "novel_allele_inference" = invisible(novel_allele_inference_project(outdir)),
            "genotype" = invisible(genotype_project(outdir)),
            "tigger_bayesian_genotype" = invisible(tigger_bayesian_genotype_project(outdir)),
            "rabhit_haplotype" = invisible(rabhit_haplotype_project(outdir)),
            "reassign_alleles" = invisible(reassign_alleles_project(outdir))
    )
    
    if (!is.null(report_params[['logo']])) {
        target_logo <- normalizePath(file.path(report_params[['outdir']],"assets", "logo.png"))
        file.copy(report_params[['logo']],
                  target_logo,
                  recursive = T, overwrite=T)
    }
    
    report_params[['outdir']] <- file.path(report_params[['outdir']],"enchantr")
    # render
    xfun::in_dir(
        outdir,
        book <- render_book(
            input="index.Rmd",
            params=report_params)
    )
    book
}

# Path to the index.Rmd template of report `name`
report_template <- function(name) {
    # Template folders are named <name>_project_files, except validate_input
    template_name <- switch(name,
                            "validate_input" = "input_validation",
                            name)
    rmd <- system.file("rstudio", "templates", "project",
                       paste0(template_name, "_project_files"), "index.Rmd",
                       package="enchantr")
    if (rmd == "") {
        stop("Report template not found for '", name, "'.")
    }
    rmd
}

# Print the parameters declared in the YAML header of a report template.
# The report description is taken from the `description` field and the
# parameter descriptions from each parameter's `label` field.
print_report_params <- function(name, rmd=report_template(name)) {
    report_params <- knitr::knit_params(readLines(rmd, warn=FALSE))

    format_default <- function(p) {
        # Defaults set with `!r` are evaluated in the current session,
        # so show the expression instead of a machine-specific value.
        if (!is.null(p$expr)) return(paste0("!r ", p$expr))
        if (is.character(p$value)) return(paste0('"', p$value, '"'))
        paste(format(p$value), collapse=", ")
    }
    format_label <- function(p) {
        label <- if (is.null(p$label)) "" else p$label
        # Drop the leading "`name`:" prefix, as the name is printed already
        label <- sub(paste0("^`?", p$name, "`?:?\\s*"), "", label)
        gsub("`", "", label)
    }

    width <- getOption("width", 80)
    description <- rmarkdown::yaml_front_matter(rmd)$description
    if (is.null(description) || !nzchar(trimws(description))) {
        description <- paste("Parameter description for enchantR report", name)
    }
    cat(strwrap(description, width=width), sep="\n")
    cat("\n")
    for (p in report_params) {
        cat(p$name, "  [default: ", format_default(p), "]\n", sep="")
        label <- format_label(p)
        if (nchar(label) > 0) {
            cat(strwrap(label, width=width - 4, indent=4, exdent=4), sep="\n")
        }
        cat("\n")
    }
    invisible(report_params)
}

#' @export
render_book <- function(input,...) {
    bookdown::render_book(
    input,
    output_file=I("index.html"),
    output_format='enchantr::immcantation',
    config_file = "_bookdown.yml",
    clean=FALSE,
    new_session=FALSE,...)
}
