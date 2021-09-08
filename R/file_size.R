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
    input <- log_table %>%
        filter(field %in% c("FILE", "ALIGNER_FILE")) %>%
        select(value) %>%
        as.character()
    if (length(input)>1) {
        stop("TODO: current version of the code expects one input file.")
    }
    output <- log_table %>%
        filter(grepl("OUTPUT",field)) %>%
        select(value) %>%
        as.character()
    if (length(output)>1) {
        stop("TODO: current version of the code expects one output file.")
    }
    
    input_size <- log_table %>%
        filter(grepl("PASS|FAIL",field)) %>%
        pull(value) %>%
        as.numeric() %>%
        sum()
    
    output_size <- log_table %>%
        filter(grepl("PASS",field)) %>%
        pull(value) %>%
        as.numeric()
    
    
    data.frame(
        "input"=input,
        "output"=output,
        "input_size"=input_size,
        "output_size"=output_size,
        "task"=task
    )
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
    log_files <-  strsplit(log_files,"[ ,]")[[1]]
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
