#' Create an Immcantation input validation project
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Input Validation
#' to create the skeleton of an Immcantation Input Validation project
#' @param  path path to the directory where the project will be created
file_size_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "file_size_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 


# consoleLogNetwork(log_file)
formatConsoleLog <- function(log_file){
    log_table <- loadConsoleLog(log_file)
    task <- unique(log_table[['task']])
    
    if (length(task)>1) {
        stop("TODO: current version of the code expects one task in each log file.")
    }
    input_df <- log_table %>%
        filter(field %in% c("FILE", "ALIGNER_FILE")) %>%
        select(value) %>%
        rename(input=value)

    if (nrow(input_df)>1) {
        stop("TODO: current version of the code expects one input file.")
    }
    
    output_df <- log_table %>%
        filter(grepl("OUTPUT",field)) %>%
        select(value) %>%
        rename(output=value)
    
    if (any(grepl("PASS|FAIL",log_table[['field']]))) {
        input_size <- log_table %>%
            filter(grepl("PASS|FAIL",field)) %>%
            summarize(input_size = sum(as.numeric(value))) %>%
            mutate(input=input_df[['input']])
    } else if ((any(grepl("RECORDS",log_table[['field']])))) {
        input_size <- log_table %>%
            filter(field == "RECORDS") %>%
            select(value) %>%
            rename(input_size=value) %>%
            mutate(input=input_df[['input']],
                   input_size=as.numeric(input_size))
    } else {
        input_size <- NA
        stop("Unknown input size")
    }

    if (any(grepl("PASS|FAIL",log_table[['field']]))) {
        if (nrow(output_df) == 1) {
            output_size <- log_table %>%
                filter(grepl("PASS",field)) %>%
                summarize(output_size=as.numeric(value)) %>%
                mutate(output=output_df[['output']])
        } else {
            stop("TODO")
        }
    } else {
        if (any(grepl("OUTPUT|SIZE",log_table[['field']])) ) {
            output_size <- log_table %>%
                select(field, value) %>%
                filter(grepl("OUTPUT|SIZE", field)) %>%
                mutate(field=gsub("([1-9]+$)","-\\1",field)) %>%
                tidyr::separate(col=field, c("type", "id"), sep="-") %>%
                tidyr::pivot_wider(id_cols=id, names_from = type) %>%
                rename(output=OUTPUT,
                       output_size=SIZE) %>%
                mutate(output_size=as.numeric(output_size)) %>%
                select(-id)
        } else {
            output_size <- NA
            stop("Unknown output size")   
        }
    }
    data.frame(
        "input"=input_df,
        "output"=output_df,
        "task"=task
    ) %>%
    left_join(input_size, by="input") %>%
    left_join(output_size, by="output")
}

getConsoleLogsGraph <- function(logs) {
    vertex_size <- bind_rows(
        logs %>% select(input,input_size) %>% rename(vertex=input, num_seqs=input_size),
        logs %>% select(output,output_size)  %>% rename(vertex=output, num_seqs=output_size),        
    ) %>%
        distinct()
    g <- graph_from_data_frame(logs[,c("input","output", "task")], directed = TRUE, vertices = vertex_size)
    g
}


plotLog <- function(g) {
    ggraph(g) + 
        geom_node_point(aes(size = num_seqs), colour = "black") +
        geom_edge_link(aes(start_cap = label_rect(node1.name), 
                           end_cap = label_rect(node2.name),
                           colour = factor(task)), 
                       arrow = arrow(type = "closed", length = unit(3, 'mm'))) 
}

# log_files <- list.files(system.file("extdata", package="enchantr"), full.names=T)
#' @export
getConsoleLogs <- function(log_files, format=c("df", "graph")) {
    format <- match.arg(format)
    #log_files <-  strsplit(log_files,"[ ,]")[[1]]
    logs <- bind_rows(lapply(log_files, formatConsoleLog))
    if (format=="df") {
        logs
    } else {
        getConsoleLogsGraph(logs)
    }
}

#' @export
plotConsoleLogs <- function(log_files) {
    logs <- getConsoleLogs(log_files, format="graph")
    lapply(decompose(logs), plotLog)
}
