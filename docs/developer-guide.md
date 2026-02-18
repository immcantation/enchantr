# Developer Guide

This document explains how the `enchantr` package is organized, how it works, how to add new reports, and how to modify the styling.

## Report Generation Workflow

The diagram below provides a high-level overview of how reports are generated in enchantr.

```mermaid
flowchart LR
    A[enchantr_report] --> B[Project Setup]
    B --> C[Template Copy]
    C --> D[Render Report]
    D --> E[HTML Output]
    
    style A fill:#e1f5ff
    style E fill:#d4edda
```

**User calls `enchantr_report()`**: The main entry point for generating reports
   ```r
   enchantr_report(name = "single_cell_qc", 
                   report_params = list(outdir = "results/"))
   ```

**Project creation**: Based on the report `name`, the corresponding `*_project()` function is called
   - Example: `single_cell_qc_project(outdir)` for single cell QC reports
   - This function copies skeleton files from `inst/rstudio/templates/project/<report_name>_project_files/` to the output directory
   - If `report_params$logo` is provided, the custom logo is copied to replace the default logo

**Report rendering**: The report is rendered using `bookdown::render_book()`
   - Input file: `index.Rmd` in the project directory
   - Output format: `enchantr::immcantation` (custom HTML format)
   - Configuration: `_bookdown.yml`


## Package Organization

```
enchantr/
├── R/                          # R source code
│   ├── report.R                # Main report generation functions
│   ├── io.R                    # Input/output utilities
│   ├── format.R                # Custom output formats
│   ├── styles.R                # ggplot2 themes
│   ├── report-1.R              # report 1 functions
│   └── ...                     # Other reports
│
├── inst/                       # Templates
│   ├── assets/                 # Static assets
│   │   ├── logo.png            # Default logo
│   │   └── style.css           # CSS stylesheet
│   ├── rstudio/templates/project/  # Report templates
│   │   ├── report-1.dcf        # Report 1 project description file
│   │   ├── report-1_project_files/
│   │   │   ├── index.Rmd           # Report main R Markdown file
│   │   │   ├── _bookdown.yml       # Bookdown config
│   │   │   ├──references.bib      # Bibliography
│   │   │   └── ...
│   │   └── ...                     # Other report templates
│   ├── scripts/                # Additional scripts
│   └── js/                     # JavaScript files
│
├── man/                        # Documentation (auto-generated)
│   └── *.Rd                    # R documentation files
│
├── docs/                       # MkDocs documentation
│   ├── topics/                 # Documentation pages
│   ├── build.R                 # Documentation build script
│   └── mkdocs.yml              # MkDocs configuration
│
├── vignettes/                  # Vignettes documentation
│   └── *.Rmd                   # Vignette files
│
└── tests/                      # Test files
    └── testthat/               # Unit tests
```

The key directories in the `enchantr` package are:

- **`R/`**: Contains R source code
  - Core functionality files (e.g., `report.R`, `io.R`, `format.R`, `styles.R`) shared across reports
  - Individual report specific files (e.g., `single_cell_qc.R`, `clonal_assignment.R`, `contamination.R`). Each of 
    this files has a `*_project` function that is used to initializes a report project using the template 
    project, and it can also contain report specific functions.
  
- **`inst/`**: Contains `.Rmd` templates and accessory files
  - `inst/assets/`: Static assets like `logo.png` and `style.css`
  - `inst/rstudio/templates/project/`: Project templates for each report type
     - one DCF file for each report. These are used to create R projects for each report type.
     - one subfolder for each report with name like `_project_files`
  - `inst/scripts/`: Additional scripts
  - `inst/js/`: JavaScript files
  
- **`man/`**: Auto-generated R documentation files (`.Rd` files)

- **`docs/`**: Documentation website source files (MkDocs format)

- **`vignettes/`**: Long-form documentation in R Markdown format


### Report Templates

Report templates are located under `inst/rstudio/templates/project/`. Each report template directory is named like `report-name_project_files`
and contains:

- **`index.Rmd`**: Main R Markdown file with:
  - YAML front matter defining:
    - Report metadata (title, date)
    - Parameters (inputs, output directory, options)
  - R code chunks for analysis and visualization
  
- **`_bookdown.yml`**: Bookdown configuration file specifying:
  - Output directory
  - Chapter/figure/table labels
  - Files to include in the book
  
- **`references.bib`**: Bibliography file for citations (if needed)

- **`references.Rmd`** or **`versions.Rmd`**: Additional chapters showing the R session information, to track versions of the packages used during the analysis.

- **Sample data files** (optional): Example input files for testing. Best if input files can be loaded from a url to reduce the size of the installed package.

### Custom Output Format

The `enchantr::immcantation` output format is defined in the `R/format.R` file and extends bookdown's HTML format with:
- Custom CSS styling (`inst/assets/style.css`)
- Custom logo
- Specific layout and theme options

## Adding a New Report

To add a new report type to enchantr, follow these steps:

### 1. Create an R file with report specific functions

Create a new file `R/<report_name>.R` with a project function:

```r
#' Create an Immcantation <Report Name> Project
#'
#' @param  path path to the directory where the project will be created
#' @export
<report_name>_project <- function(path, ...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), 
                              "rstudio", "templates", "project",
                              "<report_name>_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }

    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}
```

Add any specific analysis functions needed for the report in the same file.

### 2. Create the DCF File

Create `inst/rstudio/templates/project/<report_name>.dcf`:

```
Binding: <report_name>_project
Title: Immcantation <Report Title>
OpenFiles: index.Rmd, _bookdown.yml
```

### 3. Create the Project Files Directory

Create `inst/rstudio/templates/project/<report_name>_project_files/` with:

#### a. Example `index.Rmd` setup

```rmd
--- 
title: "Immcantation - enchantR"
subtitle: "<Report Title>"
author: "`r params$author`"
date: "Updated: `r date()`"
knit: enchantr::render_book
site: bookdown::bookdown_site
documentclass: book
bibliography: "references.bib"
biblio-style: apalike
link-citations: yes
description: "enchantR <Report Name> Report"
output: enchantr::immcantation
params:
   author:
       value: "Authors: Your Name"
   # Add report-specific parameters here
   input: 
      label: "`input`: Path to input file"
      input: file
      value: "input.tsv"
   outdir:
      label: '`outdir`: Output directory'
      input: text
      value: !r file.path(getwd(),'enchantr')
   logo:
      label: "`logo`: Path to report logo"
      input: file
      value: !r file.path("assets", "logo.png")
   echo: 
      label: '`echo`: Show code in the report.'
      input: checkbox
      value: TRUE
   cache:
      label: '`cache`: Use cached results'
      input: checkbox
      value: FALSE
---

```{r global-options, include=FALSE, cache=FALSE}
knitr::opts_chunk$set(fig.width = 7, fig.height = 4, fig.path = "figures/",
                      echo = params$echo, cache = params$cache,
                      warning = FALSE, message = FALSE,
                      eval.after = "fig.cap", out_dir = params$outdir)

# Load libraries
suppressPackageStartupMessages(library("enchantr"))
# Add other required libraries

if (!dir.exists(params[["outdir"]])) {
  dir.create(params[["outdir"]], recursive = TRUE)
}

file.copy(params$logo,
          file.path(params$outdir, "assets", "logo.png"),
          recursive = TRUE, overwrite = TRUE)
```

# Section 1

Add your analysis sections here...

```{r}
# Your analysis code
```

# References

```

#### b. Example `_bookdown.yml`

```yaml
output_dir: enchantr
delete_merged_file: true
language:
   label:
      fig: 'Figure '
      tab: 'Table '
   ui:
      chapter_name: ''
rmd_files: ["index.Rmd"]
```

#### c. Example `references.bib`

Create an empty or populated bibliography file.

#### d. Additional files

Add any required data files, additional Rmd files, or other resources.

### 4. Update `enchantr_report()`

Edit `R/report.R` to add the new report to both:

1. The function signature (allowed report names):
```r
enchantr_report <- function(name=c("validate_input", 
                                   "file_size",
                                   # ... existing reports ...
                                   "<report_name>"), 
                            report_params=list()) {
```

2. The switch statement:
```r
switch (name,
        "validate_input" = invisible(validate_input_project(outdir)),
        # ... existing cases ...
        "<report_name>" = invisible(<report_name>_project(outdir))
)
```

### 6. Test the New Report

If your report has default input data and parameters, you can simply use:

```r
# Load the package
devtools::load_all()

# Test the report
enchantr_report(name = "<report_name>", 
                report_params = list(outdir = "test_output/"))
```

### 7. Add Documentation

- Create a documentation page in `docs/topics/<report_name>_project.md`
- Update `mkdocs.yml` to include the new page in the navigation

## Modifying the CSS Template

The CSS template controls the visual styling of all enchantr reports. The main CSS file is located at `inst/assets/style.css`.

### Custom Logos

To use a custom logo without modifying CSS:

```r
enchantr_report(
  name = "single_cell_qc",
  report_params = list(
    outdir = "output/",
    logo = "path/to/custom_logo.png",
    logolink = "https://your-institution.org"
  )
)
```

The logo image dimensions should ideally be around 200-300 pixels wide for best display.

## ggplot2 Theme

In addition to the CSS for HTML reports, enchantr provides a custom ggplot2 theme for consistent plot styling. The theme is defined in `R/styles.R`. You can modify this function to change the default appearance of all plots in enchantr reports.

```r
library(ggplot2)
library(enchantr)

ggplot(data, aes(x, y)) +
  geom_point() +
  theme_enchantr()
```

## Development Workflow

Follow the conventions in [CONTRIBUTING.md](https://github.com/immcantation/immcantation/tree/master/CONTRIBUTING.md).

Key points:
- Submit your changes via pull request
- Use informative variable and function names
- Add comments for complex logic
- Follow existing code formatting
- Include roxygen2 documentation for all exported functions
- Add examples to function documentation where applicable

## Getting Help

If you have questions or need assistance:

1. Check existing [GitHub Issues](https://github.com/immcantation/enchantr/issues).
2. If your issue has not been adressed before in the issue tracker:
  a. open an new issue
  b. or send an email to the [Immcantation Group](mailto:immcantation@googlegroups.com)
