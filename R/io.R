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
            file <- plot_name
        }
        p_path <- file.path(outdir, paste0(file,".RData"))
        p_list[[plot_name]][['enchantr']][['html_caption']] <- paste0(p$labels$caption," <a href='",p_path,"'>ggplot file: ",basename(p_path),"</a>")
        save(list=names(p_list), file=p_path)
    } 
    p_list[[plot_name]]
}


#' @export
eetable <- function(df, caption=NULL) {
    element_id <- deparse1(substitute(df))
    tag <- gsub("_","-",paste0("(\\#tab:",element_id,"-table)"))
    if (!is.null(caption)) {
        caption <- paste0(tag," ",caption)   
    }
    dt <- DT::datatable(df,
                        filter="top", elementId = element_id, 
                        rownames = FALSE, fillContainer = F, 
                        options = list(scrollX = TRUE),
                        caption = htmltools::tags$caption(
                            style = 'caption-side: top; text-align: left;',
                            caption
                        ))   
    # TODO: fix caption and numbering https://stackoverflow.com/questions/49819892/cross-referencing-dtdatatable-in-bookdown
    dt
}

# In results='asis' chunk
# https://github.com/rstudio/bookdown/issues/313
# Tab. \@ref(tab:threshold-summary-table)
#' @export
print_table_caption <- function(tag, text) {
    tag <- gsub("_","-",tag)
    cat("<table>", paste0("<caption>",
                          "(#tab:",tag,"-table)",
                          " ",
                          text,
                          "</caption>"),
        "</table>", sep ="\n")
}

# db <- data.frame(
#     'id'=c(1,2),
#     'subject_id'=c("A","A")
# )
#' @export
makeLabel <- function(db, fields='id'){
    db_fields <- unique(db[,fields, drop=F])
    labels <- apply(db_fields, 2, function(label_data) {
     paste(unique(label_data),collapse="-")
    })
    paste(names(labels), labels, sep="_", collapse="_")
}

#' @export
printParams <- function(p) {
    input <- stack(p) %>%
        select(ind, values) %>%
        rename( parameter = ind,
                value = values)
    eetable(input)
}