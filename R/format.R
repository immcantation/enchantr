# https://bookdown.org/yihui/rmarkdown/format-derive.html
# https://github.com/rstudio/bookdown/blob/main/R/gitbook.R
# https://github.com/rstudio/bookdown/tree/main/inst/resources/gitbook
# https://github.com/rstudio/bookdown/blob/main/inst/templates/gitbook.html
# https://github.com/rstudio/bookdown/blob/main/inst/rstudio/templates/project/resources/gitbook/_output.yml
# https://github.com/rstudio/bookdown/issues/61#issuecomment-200996786

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
    if (!is.null(user_args)) {
      immcantation_config <- modifyList(
        immcantation_config, 
        user_args, 
        keep.null = T) 
    }
    
    # overwrite bookdown:::gitbook_dependency
    assignInNamespace("gitbook_dependency",gitbook_dependency,ns="bookdown")
    
    # call the base html_document function
    b <- bookdown::gitbook(
                      toc_depth= immcantation_config[['config']][['toc']][['depth']],
                      css="assets/style.css",
                      fig_caption=immcantation_config[['config']][['fig_caption']],
                      highlight = immcantation_config[['highlight']],
                      keep_md = immcantation_config[['config']][['keep_md']],
                      config=c(immcantation_config[['config']],
                               immcantation_config['fontsettings']
                               ),
                      code_folding="hide", self_contained=TRUE, split_by = "none" 
                      )
    invisible(b)
}


# This function overwrites bookdown:::gitbook_dependency
# to use a local copy of fuse.js.
# The original function has a hardcoded url to download fuse. fuse is the search
# engine for the reports. This is a proble for users in controlled environments 
# without internet access because the can't render the reports. The error is:
# "Could not fetch https://cdn.jsdelivr.net/npm/fuse.js@6.4.6/dist/fuse.min.js"
gitbook_dependency = function(table_css, config = list()) {
    assets = bookdown:::bookdown_file('resources', 'gitbook')
    owd = setwd(assets); on.exit(setwd(owd), add = TRUE)
    app = if (file.exists('js/app.min.js')) 'app.min.js' else 'app.js'
    # TODO: download and a local copy of fuse.js?
    fuse = htmltools::htmlDependency(
        'fuse', '6.4.6', c('file' = system.file("js","fuse.js", package = "enchantr", lib.loc = NULL,
                                                mustWork = FALSE)),
        script = 'fuse.js@6.4.6'
    )
    if (is.logical(config$search)) {
        lunr = FALSE
        if (!config$search) fuse = NULL
    } else {
        # use fuse as the search engine by default
        lunr = identical(config$search$engine, 'lunr')
    }
    list(bookdown:::jquery_dependency(), fuse, htmltools::htmlDependency(
        'gitbook', '2.6.7', src = assets,
        stylesheet = file.path('css', c(
            'style.css', if (table_css) 'plugin-table.css', 'plugin-bookdown.css',
            'plugin-highlight.css', 'plugin-search.css', 'plugin-fontsettings.css',
            'plugin-clipboard.css'
        )),
        script = file.path('js', c(
            app, if (lunr) 'lunr.js', 'clipboard.min.js', 'plugin-search.js', 'plugin-sharing.js',
            'plugin-fontsettings.js', 'plugin-bookdown.js', 'jquery.highlight.js',
            'plugin-clipboard.js'
        ))
    ))
}
