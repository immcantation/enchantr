# Report names kept for backwards compatibility. Each resolves to the unified
# "genotype" report with its `method` parameter preset.
.GENOTYPE_REPORT_ALIASES <- c(
    piglet_genotype = "allele_based"
)

#' Render an Immcantation project
#'
#' @param  name report name
#' @param  report_params params list, needed by report. must include outdir
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
                                   "piglet_genotype",
                                   "reassign_alleles"), report_params=list()) {

    name <- match.arg(name)

    # Resolve a deprecated genotype report name to the unified report.
    if (name %in% names(.GENOTYPE_REPORT_ALIASES)) {
        method <- .GENOTYPE_REPORT_ALIASES[[name]]
        message("Report '", name, "' is deprecated; use enchantr_report('genotype', ",
                "list(method = '", method, "', ...)) instead.")
        if (is.null(report_params[['method']])) {
            report_params[['method']] <- method
        }
        name <- "genotype"
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