# https://bookdown.org/yihui/rmarkdown/format-derive.html
# https://github.com/rstudio/bookdown/blob/main/R/gitbook.R
# https://github.com/rstudio/bookdown/tree/main/inst/resources/gitbook
# https://github.com/rstudio/bookdown/blob/main/inst/templates/gitbook.html
# https://github.com/rstudio/bookdown/blob/main/inst/rstudio/templates/project/resources/gitbook/_output.yml

#' @importFrom rmarkdown output_format knitr_options pandoc_options
#' @export
immcantation <- function(...) {
    
    user_args <- list(...)
    
    # copy css, don't overwrite user provided logo.png
    file.copy(
        system.file("assets", package = "enchantr"),
        getwd(),
        recursive = T, overwrite = F)
    
    enchantr_dir <- file.path(getwd(), "enchantr")
    if (!dir.exists(enchantr_dir)) {
        dir.create(enchantr_dir, recursive = T)
    }
    
    file.copy(
        "assets",
        enchantr_dir,
        recursive = T, overwrite = F)

    immcantation_config <- list(
        "config" = list (
          "fig_caption" = T,
          "toc" = list(
            "scroll_highlight" = T,
            "collapse" = "subsection",
            "depth"=4,
            "before"="<li class='toc-logo'><a href='http://immcantation.readthedocs.io' target='blank'><img src='./assets/logo.png'></a></li>"
          ),
          "download"=F,
          "sharing"=F,
          "keep_md"=F
        ),
        "fontsettings"=list(
          "theme"="white",
          "familiy"="sans",
          "size"=1
        ),
        "highlight"="pygments"
    )
    if (!is.null(user_args)){
      immcantation_config <- modifyList(
        immcantation_config, 
        user_args, 
        keep.null = T) 
    }
    
    # call the base html_document function
    b <- bookdown::gitbook(
                      toc_depth= immcantation_config[['config']][['toc']][['depth']],
                      css="assets/style.css",
                      fig_caption=immcantation_config[['config']][['fig_caption']],
                      highlight = immcantation_config[['highlight']],
                      keep_md = immcantation_config[['config']][['keep_md']],
                      config=c(immcantation_config[['config']],immcantation_config['fontsettings']),
                      code_folding="show"
                      )
    invisible(b)
}