#' @export
theme_enchantr <- function() {
    ggplot2::theme_bw(base_family = "ArialMT") +
        theme(strip.background = element_blank(),
              plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(family = "ArialMT")
              )
}
