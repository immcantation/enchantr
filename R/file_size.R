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
        stop("more than one task in the same log")
    }
    
    if (task == "ParseDb-split") {
        file_field <- grepl("FILE",log_table[['field']])
        if (sum(file_field)>1) {
            stop("Can't deal with multiple input files.")
        }
        log_table[['field']][file_field] <- 'FILE1'
        records_field <- grepl("RECORDS",log_table[['field']])
        log_table[['field']][records_field] <- 'RECORDS1'
        parts <- as.numeric(log_table[['value']][log_table[['field']]=="PARTS"])
        i <- 2
        while (i <= parts) {
            log_table <- bind_rows(log_table,
                      log_table[which(file_field),]   
            )
            log_table[['field']][nrow(log_table)] <- paste0("FILE",i)
            i <- i + 1
        }
        
        
    }
    
    if (task == "MakeDB-igblast") {
        log_table <- log_table %>%
            filter(field != "SEQ_FILE")
    }
    
    .getFileID <- function(x) {
        sapply(x, function(i){
            if (grepl("[0-9]$",i)){
                gsub("(.*)([0-9]+)$","\\2",i)  
            } else {
                0
            }
        })
        
    }
    log_table <- log_table %>%
        filter(    grepl("ALIGNER_FILE",field) |
                       grepl("SEQ_FILE",field) |
                       grepl("FILE[0-9]*",field) |
                       grepl("PASS",field) |
                       grepl("FAIL",field) |
                       grepl("RECORDS",field) |
                       grepl("SIZE",field) |
                       grepl("OUTPUT[0-9]*",field)                     
        ) %>%
        mutate(file_id=.getFileID(field)) %>%
        mutate(
            field_type = case_when(grepl("OUTPUT", field) ~ "output", 
                                   grepl("FILE", field) ~ "input",
                                   grepl("ALIGNER_FILE", field) ~ "input",
                                   grepl("SEQ_FILE", field) ~ "input",
                                   grepl("PASS", field) ~ "output_size",
                                   grepl("RECORDS", field) ~ "input_size",
                                   grepl("SIZE", field) ~ "output_size",
                                   grepl("FAIL", field) ~ "fail")
        ) %>%
        select(-field, -step) 
    
    log_table_input <- log_table %>%
        filter(grepl("input", field_type))
    input_fields <- log_table_input[["field_type"]] == "input"
    
    if (!all(input_fields)) {
        log_table_input <- log_table_input[input_fields,] %>%
            left_join(
                log_table_input <- log_table_input[!input_fields,] %>%
                    pivot_wider(id_cols=file_id, 
                                values_from=value, names_from=field_type)
            )
    }
    log_table_input <- log_table_input %>%
        rename(input=value)
    
    log_table_output <- log_table %>%
        filter(!grepl("input", field_type)) %>%
        pivot_wider(id_cols=file_id, 
                    values_from=value, names_from=field_type)
    
    log_table <- log_table_input %>%
        left_join(log_table_output)
    
    if (!"input_size" %in% colnames(log_table)) {
        if (all(c("output_size", "fail") %in% colnames(log_table))) {
            log_table <- log_table %>%
                group_by(file_id) %>%
                rowwise() %>%
                mutate(input_size=as.numeric(output_size)+as.numeric(fail))   
        } else {
            log_table[['input_size']] <- NA
        }
    } 
    log_table[['task']] <- task
    log_table[['input_size']] <- as.numeric(log_table[['input_size']])
    log_table[['output_size']] <- as.numeric(log_table[['output_size']])
    log_table %>%
        ungroup() %>%
        select(input, output, task, input_size, output_size) %>%
        as.data.frame(stringsAsFactors = F)
    
}

getConsoleLogsGraph <- function(logs) {
    vertex_size <- bind_rows(
        logs %>% select(input,input_size) %>% rename(vertex=input, num_seqs=input_size),
        logs %>% select(output,output_size)  %>% rename(vertex=output, num_seqs=output_size),        
    ) %>%
    distinct() 
    # rm duplicate vertex names
    dup_v <- duplicated(vertex_size$vertex)
    if (any(dup_v)) {
        rm_me <- vertex_size$vertex %in% vertex_size$vertex[duplicated(vertex_size$vertex)] &
            is.na(vertex_size$num_seqs)
        vertex_size <- vertex_size[!rm_me,,drop=F]
    }
    g <- graph_from_data_frame(logs[,c("input","output", "task")], 
                               directed = TRUE, 
                               vertices = vertex_size)
    g
}


plotLog <- function(g) {
    ggraph(g) + 
        #geom_node_point(aes(size = num_seqs), colour = "black") +
        geom_node_label(aes(label=paste0(name,": ", num_seqs, " sequences")), 
                        colour = "black", label.size=0) +
        geom_edge_link(aes(start_cap = label_rect(node1.name), 
                           end_cap = label_rect(node2.name),
                           # colour = factor(task),
                           label=task), 
                       color="grey50",
                       angle_calc = "across",
                       arrow = arrow(type = "closed", 
                                     length = unit(3, 'mm'))) +
        theme_enchantr() +
        theme(
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.border = element_blank()
        )
}

plotWorkflow <- function(g) {
    root_nodes <- which(degree(g, mode = "in") == 0)
    g <- add_vertices(g, 1, attr = list("name"="Input"))
    input_node <- vcount(g)
    for (i in root_nodes) {
        g <- add_edges(g, c(input_node,i), task="Input")
    }
    g2 <- make_line_graph(g) 
    V(g2)$name <- E(g)$task
    g3 <- graph_from_edgelist(unique(as_edgelist(g2)))
    p <- ggraph(g3, layout = "auto") +
    # geom_edge_link0(aes(col=1,edge_width=2), arrow = arrow(length = unit(8, 'mm')))+
        geom_edge_link(color="darkblue",edge_width=3, arrow = arrow(length = unit(6, 'mm')))+
        geom_node_point(shape=21,col="white",fill="black",size=5,stroke=1)+
        geom_node_label(aes(label = name), repel = FALSE) +
        theme_graph(plot_margin = margin(5, 50, 5, 40)) +
        coord_cartesian(clip = "off") +
        theme_enchantr() +
        theme(
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.border = element_blank()
        )
    p
}

# log_files <- list.files(system.file("extdata", package="enchantr"), full.names=T)
#' @export
getConsoleLogs <- function(log_files, format=c("df", "graph")) {
    format <- match.arg(format)
    #log_files <-  strsplit(log_files,"[ ,]")[[1]]
    logs <- bind_rows(lapply(log_files, enchantr:::formatConsoleLog))
    if (format=="df") {
        logs
    } else {
        getConsoleLogsGraph(logs)
    }
}

#' @export
plotConsoleLogs <- function(log_files, style=c("decompose", "workflow")) {
    style <- match.arg(style)
    logs <- getConsoleLogs(log_files, format="graph")
    if (style == "decompose") {
        lapply(decompose(logs), enchantr:::plotLog)
    } else {
        plotWorkflow(logs)
    }

}
