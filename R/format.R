# https://bookdown.org/yihui/rmarkdown/format-derive.html
# https://github.com/rstudio/bookdown/blob/main/R/gitbook.R
# https://github.com/rstudio/bookdown/tree/main/inst/resources/gitbook
# https://github.com/rstudio/bookdown/blob/main/inst/templates/gitbook.html
# https://github.com/rstudio/bookdown/blob/main/inst/rstudio/templates/project/resources/gitbook/_output.yml

#' @importFrom rmarkdown output_format knitr_options pandoc_options
#' @export
immcantation = function(...) {
    
    file.copy(
        system.file("assets", package = "enchantr"),
        getwd(),
        recursive = T)
    # TODO: there must be a better way to do this
    file.copy(
        system.file("assets", package = "enchantr"),
        file.path(getwd(), "enchantr"),
        recursive = T)
    user_args <- list(...)
    save(user_args, file="user_args.RData")
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
    save(immcantation_config, file="immcantation_config_default.RData")
    if (!is.null(user_args)){
      immcantation_config <- modifyList(
        immcantation_config, 
        user_args, 
        keep.null = T) 
    }
    
    save(immcantation_config, file="immcantation_config.RData")
    # call the base html_document function
    b <- bookdown::gitbook(
                      toc_depth= immcantation_config[['config']][['toc']][['depth']],
                      css="assets/style.css",
                      fig_caption=immcantation_config[['config']][['fig_caption']],
                      highlight = immcantation_config[['highlight']],
                      keep_md = immcantation_config[['config']][['keep_md']],
                      config=c(immcantation_config[['config']],immcantation_config['fontsettings'])
                      )
    # file.copy(
    #   file.path("assets", "logo.png"),
    #   file.path(getwd(), "enchantr", "assets"),
    #   recursive = T)
    save(b,file="b.RData") 
    invisible(b)
}