#### enchantr_report(help=TRUE) ####

report_names <- eval(formals(enchantr_report)$name)

# Write a minimal report template with the given YAML header lines
write_template <- function(header) {
    rmd <- tempfile(fileext=".Rmd")
    writeLines(c("---", header, "---", "", "Report body"), rmd)
    rmd
}

test_that("help prints the report description and parameters", {
    out <- capture.output(
        res <- enchantr_report("collapse_duplicates", help=TRUE)
    )
    expect_equal(out[1], "enchantr collapse duplicates report")
    expect_true(any(grepl('^input  \\[default: "input.tsv"\\]$', out)))
    expect_true(any(grepl('^collapseby  \\[default: "sample_id"\\]$', out)))
    expect_true(any(grepl("^nproc  \\[default: 1\\]$", out)))
    expect_true(any(grepl("^echo  \\[default: TRUE\\]$", out)))
    # Returns the parsed parameters
    expect_true(all(c("input", "collapseby", "num_fields", "text_fields",
                      "outdir") %in% names(res)))
})

test_that("help returns the parameters invisibly", {
    capture.output(vis <- withVisible(
        enchantr_report("collapse_duplicates", help=TRUE)
    ))
    expect_false(vis$visible)
})

test_that("help does not create outdir or render the report", {
    outdir <- file.path(tempdir(), "enchantr_help_test")
    unlink(outdir, recursive=TRUE)
    capture.output(
        enchantr_report("collapse_duplicates",
                        report_params=list(outdir=outdir), help=TRUE)
    )
    expect_false(dir.exists(outdir))
})

test_that("help works for every report", {
    for (name in report_names) {
        expect_true(file.exists(enchantr:::report_template(name)), info=name)
        capture.output(res <- enchantr_report(name, help=TRUE))
        expect_gt(length(res), 0, label=name)
    }
})

test_that("validate_input maps to the input_validation template", {
    expect_match(enchantr:::report_template("validate_input"),
                 "input_validation_project_files/index.Rmd$")
})

test_that("unknown report names give an error", {
    expect_error(enchantr_report("not_a_report", help=TRUE))
    expect_error(enchantr:::report_template("not_a_report"),
                 "Report template not found")
})

test_that("missing description prints a default header", {
    rmd <- write_template(c(
        "params:",
        "   nproc:",
        "      label: \"`nproc`: Number of cores\"",
        "      value: 1"
    ))
    on.exit(unlink(rmd))
    out <- capture.output(enchantr:::print_report_params("my_report", rmd=rmd))
    expect_equal(out[1], "Parameter description for enchantR report my_report")
})

test_that("empty description prints a default header", {
    rmd <- write_template(c(
        "description: \"\"",
        "params:",
        "   nproc:",
        "      value: 1"
    ))
    on.exit(unlink(rmd))
    out <- capture.output(enchantr:::print_report_params("my_report", rmd=rmd))
    expect_equal(out[1], "Parameter description for enchantR report my_report")
})

test_that("labels and defaults are formatted", {
    rmd <- write_template(c(
        "description: \"My test report\"",
        "params:",
        "   input:",
        "      label: \"`input`: Path to `tsv` file\"",
        "      value: \"input.tsv\"",
        "   outdir:",
        "      label: \"`outdir`: Output directory\"",
        "      value: !r file.path(getwd(), 'enchantr')",
        "   nolabel:",
        "      value: 3",
        "   plainlabel:",
        "      label: \"plainlabel: No backticks\"",
        "      value: FALSE"
    ))
    on.exit(unlink(rmd))
    out <- capture.output(enchantr:::print_report_params("my_report", rmd=rmd))
    expect_equal(out[1], "My test report")
    # Name prefix and backticks are removed from labels
    expect_true('input  [default: "input.tsv"]' %in% out)
    expect_true("    Path to tsv file" %in% out)
    # `!r` defaults show the expression, not the evaluated value
    expect_true("outdir  [default: !r file.path(getwd(), 'enchantr')]" %in% out)
    expect_false(any(grepl(getwd(), out, fixed=TRUE)))
    # Parameters without label print only name and default
    i <- which(out == "nolabel  [default: 3]")
    expect_length(i, 1)
    expect_equal(out[i + 1], "")
    expect_true("plainlabel  [default: FALSE]" %in% out)
    expect_true("    No backticks" %in% out)
})

test_that("long descriptions are wrapped to the console width", {
    rmd <- write_template(c(
        "params:",
        "   input:",
        paste0("      label: \"`input`: ", strrep("word ", 40), "\""),
        "      value: \"x\""
    ))
    on.exit(unlink(rmd))
    old <- options(width=60)
    on.exit(options(old), add=TRUE)
    out <- capture.output(enchantr:::print_report_params("my_report", rmd=rmd))
    label_lines <- grep("^    word", out, value=TRUE)
    expect_gt(length(label_lines), 1)
    expect_true(all(nchar(label_lines) <= 60))
})
