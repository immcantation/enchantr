# Developer Guide

This document explains how the `enchantr` package is organized, how it works, how to add new reports, and how to modify the styling.

[TOC]

## Package Organization

### Directory Structure

The key directories in the enchantr package are:

- **`R/`**: Contains R source code
  - Core functionality files (e.g., `report.R`, `io.R`, `format.R`, `styles.R`) shared across reports
  - Individual report specific files (e.g., `single_cell_qc.R`, `clonal_assignment.R`, `contamination.R`). Each of 
    this files has a `*_project` function that is used to initializes a report project using the template 
    project, and it can also contain report specific functions.
  
- **`inst/`**: Contains `.Rmd` templates and accessory files
  - `inst/assets/`: Static assets like `logo.png` and `style.css`
  - `inst/rstudio/templates/project/`: Project templates for each report type
     - one DCF file for each report. These are used to create R projects for each report type.
     - one subfolder for each report
  - `inst/scripts/`: Additional scripts
  - `inst/js/`: JavaScript files
  
- **`man/`**: Auto-generated R documentation files (`.Rd` files)

- **`docs/`**: Documentation website source files (MkDocs format)

- **`vignettes/`**: Long-form documentation in R Markdown format

## Report Generation Workflow

1. **User calls `enchantr_report()`**: The main entry point for generating reports
   ```r
   enchantr_report(name = "single_cell_qc", 
                   report_params = list(outdir = "results/"))
   ```

2. **Project creation**: Based on the report `name`, the corresponding `*_project()` function is called
   - Example: `single_cell_qc_project(outdir)` for single cell QC reports
   - This function copies skeleton files from `inst/rstudio/templates/project/<report_name>_project_files/` to the output directory

3. **Logo customization** (optional): If `report_params$logo` is provided, the custom logo is copied to replace the default logo

4. **Report rendering**: The report is rendered using `bookdown::render_book()`
   - Input file: `index.Rmd` in the project directory
   - Output format: `enchantr::immcantation` (custom HTML format)
   - Configuration: `_bookdown.yml`



### Report Template Structure

Each report template directory contains:

- **`index.Rmd`**: Main R Markdown file with YAML front matter defining:
  - Report metadata (title, author, date)
  - Parameters (inputs, output directory, options)
  - R code chunks for analysis and visualization
  
- **`_bookdown.yml`**: Bookdown configuration file specifying:
  - Output directory
  - Chapter/figure/table labels
  - Files to include in the book
  
- **`references.bib`**: Bibliography file for citations

- **`references.Rmd`** or **`versions.Rmd`**: Additional chapters

- **Sample data files** (optional): Example input files for testing

### Custom Output Format

The `enchantr::immcantation` output format is defined in the `R/format.R` file and extends bookdown's HTML format with:
- Custom CSS styling (`inst/assets/style.css`)
- Custom logo
- Specific layout and theme options

## Adding a New Report

To add a new report type to enchantr, follow these steps:

### 1. Create the R Function

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

#### a. `index.Rmd`

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

#### b. `_bookdown.yml`

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

#### c. `references.bib`

Create an empty or populated bibliography file.

#### d. Additional files

Add any required data files, additional Rmd files, or other resources.

### 4. Update `enchantr_report()`

Edit `R/report.R` to add the new report to:

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

### 5. Document the Function

Run `devtools::document()` to generate the `.Rd` file in `man/`.

### 6. Test the New Report

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

The CSS template controls the visual styling of all enchantr reports.

### Location

The main CSS file is located at:
```
inst/assets/style.css
```

### CSS Structure

The style.css file contains styles for:

1. **Fonts**: Google Fonts imports and font family definitions
   ```css
   @import url('https://fonts.googleapis.com/css?family=...');
   
   .book.font-family-1 {
     font-family: 'Karla', arial, sans-serif;
   }
   
   h1, h2, h3, h4 {
     font-family: 'Lora', arial, sans-serif;
   }
   ```

2. **Page Structure**: Background colors, padding, margins
   ```css
   .page-inner {
     padding-top: 10px !important;
     padding-bottom: 0px !important;
     background-color: #3fb5bd;
   }
   ```

3. **Links**: Color and styling
   ```css
   .book .book-body .page-wrapper .page-inner section.normal a {
     color: #157f8e;
   }
   ```

4. **Headers**: Size, color, and spacing

5. **Tables**: Borders, padding, and alignment

6. **Code blocks**: Background, borders, and text formatting

7. **Figures and captions**: Sizing and positioning

8. **Navigation elements**: Sidebar, top bar, buttons

### Making Changes

1. **Edit the CSS file**:
   ```bash
   # Open the file
   vim inst/assets/style.css
   ```

2. **Modify styles**: Use standard CSS syntax
   ```css
   /* Example: Change link color */
   .book .book-body .page-wrapper .page-inner section.normal a {
     color: #ff0000;  /* Red links */
   }
   
   /* Example: Change header font */
   h1, h2, h3, h4 {
     font-family: 'Arial', sans-serif;
     color: #333333;
   }
   
   /* Example: Change background color */
   .page-inner {
     background-color: #ffffff;  /* White background */
   }
   ```

3. **Test changes**:
   ```r
   # Reinstall package to include updated CSS
   devtools::install()
   
   # Generate a test report
   enchantr_report(name = "single_cell_qc", 
                   report_params = list(outdir = "test_css/"))
   
   # Open the report in a browser to view changes
   browseURL("test_css/enchantr/index.html")
   ```

### CSS Tips

- **Browser Developer Tools**: Use your browser's developer tools (F12) to inspect elements and test CSS changes interactively before editing the file

- **Specificity**: Some styles may need `!important` to override default bookdown styles
  ```css
  .element {
    property: value !important;
  }
  ```

- **Responsive Design**: Consider different screen sizes
  ```css
  @media screen and (max-width: 768px) {
    /* Styles for mobile devices */
  }
  ```

- **Color Schemes**: Maintain consistency across the report by defining color variables (if using CSS preprocessors) or keeping a color palette comment

- **Testing**: Always test reports with different content types (tables, figures, code blocks) to ensure styles work universally

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

In addition to the CSS for HTML reports, enchantr provides a custom ggplot2 theme for consistent plot styling.

### Usage

```r
library(ggplot2)
library(enchantr)

ggplot(data, aes(x, y)) +
  geom_point() +
  theme_enchantr()
```

### Customization

The theme is defined in `R/styles.R`:

```r
theme_enchantr <- function() {
    ggplot2::theme_bw(base_family = "ArialMT") +
        theme(strip.background = element_blank(),
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(family = "ArialMT")
              )
}
```

You can modify this function to change the default appearance of all plots in enchantr reports.

## Development Workflow

### Making Changes

1. **Clone the repository**:
   ```bash
   git clone https://github.com/immcantation/enchantr.git
   cd enchantr
   ```

2. **Create a new branch** for your feature:
   ```bash
   git checkout -b feature-name
   ```

3. **Make your changes** following the guidelines above

4. **Document your changes**:
   ```r
   devtools::document()
   ```

5. **Test your changes**:
   ```r
   devtools::test()
   devtools::check()
   ```

6. **Commit and push**:
   ```bash
   git add .
   git commit -m "Description of changes"
   git push origin feature-name
   ```

7. **Create a pull request** on GitHub

### Code Style

Follow the conventions in [CONTRIBUTING.md](https://github.com/immcantation/immcantation/tree/master/CONTRIBUTING.md).

Key points:
- Use informative variable and function names
- Add comments for complex logic
- Follow existing code formatting
- Include roxygen2 documentation for all exported functions
- Add examples to function documentation where applicable

## Additional Resources

- **Bookdown documentation**: https://bookdown.org/yihui/bookdown/
- **R Markdown documentation**: https://rmarkdown.rstudio.com/
- **AIRR Standards**: https://docs.airr-community.org/
- **Immcantation framework**: https://immcantation.readthedocs.io/

## Getting Help

If you have questions or need assistance:

1. Check existing [GitHub Issues](https://github.com/immcantation/enchantr/issues)
2. Contact the [Immcantation Group](mailto:immcantation@googlegroups.com)
3. Join the [Immcantation News](https://groups.google.com/g/immcantation-news) Google Group
