#' @export
theme_enchantr <- function() {
    ggplot2::theme_bw() +
        theme(strip.background=element_blank(),
              plot.background=element_blank(),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank()
              )
}

