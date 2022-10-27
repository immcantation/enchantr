#' eeplot: a convenience function to save ggplot figures
#' 
#' This function takes a \code{ggplot}, and saves it in an .RData file in
#' the \code{outdir} directory. It returns the same input \code{p} with an 
#' additional field,\code{enchantr}, with an \code{html_caption}
#' that can be used in html reports, to provide a download link to the figure file.
#' 
#' @param p ggplot figure
#' @param outdir directory where the output file will be saved.
#' @param file   filename. If null, \code{p} will be used.
#' @param caption catption for the image. Will be updated to add the path to
#'                the output file.
#' @param ...   Additional objects to be save in the same output file               
#' @examples 
#' diamonds_plot <- ggplot(ggplot2::diamonds, aes(carat)) + geom_histogram() +
#'     labs(title = "Title of the plot",
#'          subtitle = "Subtitle of the plot",
#'          caption = "This is the caption")
#' eeplot(diamonds_plot, outdir=tempdir(), file="diamonds-plot")
#' @export
eeplot <- function(p, outdir=NULL, file=NULL, caption=NULL, ... ) {
    # This hack is for the plot to maintain the original name in the file,
    # and load it with the same name, not 'p'.
    plot_name <- deparse1(substitute(p))
    # p[['enchantr']][['html_caption']] <- p$labels$caption
    p[['enchantr']][['html_caption']] <- caption
    p_list <- list("p"=p)
    names(p_list) <- plot_name
    assign(plot_name, p)
    if (!is.null(outdir)) {
        if (is.null(file)) {
            file <- plot_name
        }
        p_path <- paste0(file,".RData")
        p_list[[plot_name]][['enchantr']][['html_caption']] <- paste0(p$labels$caption," <a href='",p_path,"'>ggplot file: ",basename(p_path),"</a>")
        p_path <- file.path(outdir, paste0(file,".RData"))
        var_name <- paste(plot_name,"extra_objects", sep="_")
        assign(var_name, list(...))
        save(list=c(names(p_list),var_name), file=p_path)
    } 
    p_list[[plot_name]]
}


#' eeplot: a convenience function to save tables
#' 
#' This function takes a \code{data.frame}, and saves it as a .tsv file in
#' the \code{outdir} directory. 
#' 
#' @param df data.frame
#' @param caption caption for the table
#' @param outdir directory where the output table be saved.
#' @param file   file name (without extension). If null, will use the name of the
#'               input data.frame.
#' @param file   filename. If null, \code{p} will be used.
#' @param show_max Number of lines to show. All data will be saved in the .tsv file.
#'
#' @return It returns a list with a \code{DT::datatable} and the caption, 
#'         updated with additional text with a link to the destination file.
#' @seealso  See \link{print_table_caption}, the use of which is required for
#'           correct placement of the caption.
#' @export
eetable <- function(df, caption=NULL, outdir=NULL, file=NULL, show_max=NULL) {
    element_id <- deparse1(substitute(df))
    
    # caption and numbering https://stackoverflow.com/questions/49819892/cross-referencing-dtdatatable-in-bookdown
    # tag <- gsub("_","-",paste0("(\\#tab:",element_id,")"))
    tag <- ""
    if (!is.null(caption)) {
        caption <- paste0(tag," ",caption)   
    }
    if (!is.null(outdir)) {
        if (is.null(file)) {
            file <- element_id
        }
        #tab_path <- file.path(outdir, paste0(file,".tsv"))
        tab_path <- paste0(file,".tsv")
        caption <- paste0(caption, " File can be found here: <a href='",tab_path,"'> ",basename(tab_path),"</a>")
        write.table(df, 
                    file = file.path(outdir,tab_path), 
                    sep="\t", quote = F, row.names = F)
    } 
    if (!is.null(show_max)) {
        df <- df[1:min(show_max,nrow(df)),]
    }
    dt <- DT::datatable(df,
                        filter="top", elementId = element_id, 
                        rownames = FALSE, fillContainer = F, 
                        options = list(
                            scrollX = TRUE)#,
                        # caption = htmltools::tags$caption(
                        #     style = 'caption-side: top; text-align: left;',
                        #     caption
                        # )
    )   
    # https://stackoverflow.com/questions/70868546/r-markdown-printing-datatable-from-inside-the-function
    # https://stackoverflow.com/questions/30509866/for-loop-over-dygraph-does-not-work-in-r
    # doesn't work, prints whit space
    # print(htmltools::tagList(list(dt)))
    list(
        "table"=dt,
        "caption"=caption
    )
}

# Hack to print table captions correctly
# Must be used right after eetable, in a results='asis' chunk
# https://github.com/rstudio/bookdown/issues/313
# Tab. \@ref(tab:threshold-summary)
# No _! _ converted to -
# tag  the name of the data.frame
# text the caption
#' @export
print_table_caption <- function(tag, text) {
    # rm latex tax if present
    if (!is.null(text)) {
        if  (grepl("^\\(", text)) {
            text <- sub(".*\\)\\s+","", text)
        }   
    }
    tag <- gsub("_","-",tag)
    cat("<table>", paste0("<caption>",
                          "(#tab:",tag,")",
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
    # too long
    #paste(names(labels), labels, sep="_", collapse="_")
    paste(labels, collapse="_")
}

#' @export
printParams <- function(p) {
    input <- stack(p) %>%
        select(ind, values) %>%
        rename( parameter = ind,
                value = values)
    eetable(input)$table
}