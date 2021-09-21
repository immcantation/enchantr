#' eeplot: a convenience function to save ggplot figures
#' 
#' This function takes a \code{ggplot}, and saves it in an .RData file in
#' the \cod{outdir} directory. It returns the same input \code{p} with an 
#' additional field,\code{enchantr}, with an \code{html_caption}
#' that can be used in html reports, to provide a download link to the figure file.
#' 
#' @param p ggplot figure
#' @param outdir directory where the output file will be saved.
#' @param file   filename. If null, \code{p} will be used.
#'
#' @examples 
#' diamonds_plot <- ggplot(ggplot2::diamonds, aes(carat)) + geom_histogram() +
#'     labs(title = "Title of the plot",
#'          subtitle = "Subtitle of the plot",
#'          caption = "This is the caption")
#' eeplot(diamonds_plot, outdir=tempdir(), file="diamonds-plot")
#' @export
eeplot <- function(p, outdir=NULL, file=NULL) {
    # This hack is for the plot to maintain the original name in the file,
    # and load it with the same name, not 'p'.
    plot_name <- deparse1(substitute(p))
    p[['enchantr']][['html_caption']] <- p$labels$caption
    p_list <- list("p"=p)
    names(p_list) <- plot_name
    assign(plot_name, p)
    if (!is.null(outdir)) {
        if (is.null(file)) {
            file <- deparse1(substitute(p))
        }
        p_path <- file.path(outdir, paste0(file,".RData"))
        p_list[[plot_name]][['enchantr']][['html_caption']] <- paste0(p$labels$caption," <a href='",p_path,"'>ggplot file: ",basename(p_path),"</a>")
        save(list=names(p_list), file=p_path)
    } 
    p_list[[plot_name]]
}


#' @export
eetable <- function(df, caption) {
    element_id <- deparse1(substitute(p))
    caption <- paste0(paste0("(\\#tab:",element_id,"-table"),caption)
    DT::datatable(df,
                  filter="top", elementId = element_id, 
                  rownames = FALSE, fillContainer = F, 
                  options = list(scrollX = TRUE),
                  caption=caption)
    
}

