#' Create an Immcantation input validation project
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Input Validation
#' to create the skeleton of an Immcantation Input Validation project
#' @param  path path to the directory where the project will be created
file_size_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                              "templates", "project",
                              "file_size_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}


#' Format Console Log Table for File Size Tracking
#'
#' formatConsoleLog processes a console log file or a data frame containing log information,
#' and formats it to extract and summarize input/output file names and their corresponding
#' record and size information for downstream analysis.
#'
#' @param log_file Character string specifying the path to a log file, or a data.frame
#'   containing log information (for testing purposes).
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{input}{Input file name(s)}
#'     \item{output}{Output file name(s)}
#'     \item{task}{Task name}
#'     \item{input_size}{Number of input records (numeric)}
#'     \item{output_size}{Number of output records (numeric)}
#'   }
#'
#' @details
#' The function supports logs from different tasks (e.g., "ParseDb-split", "MakeDB-igblast"),
#' handles multiple input/output files, and computes record counts for each file.
#' It also normalizes field names and merges related log entries for easier downstream use.
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' # Using a log file path
#' formatConsoleLog("path/to/log_file.log")
#'
#' # Using a data.frame
#' log_df <- structure(
#'    list(step = c(1, 1, 1, 1, 1, 1, 1, 1, 1), 
#'          task = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
#'                           levels = "ParseDb-split", class = "factor"), 
#'          field = c("FILE", "FIELD", "NUM_SPLIT", "OUTPUT1", "OUTPUT2", 
#'                    "RECORDS", "PARTS", "SIZE1", "SIZE2"), 
#'          value = c("Sample8_quality-pass.tsv", "productive", "None", 
#'                    "Sample8_productive-F.tsv", "Sample8_productive-T.tsv", 
#'                    "92", "2", "3", "89")), 
#'     class = "data.frame", 
#'     row.names = 3:11)
#' formatConsoleLog(log_df)
#' }
#'
formatConsoleLog <- function(log_file){
    
    # allow passing either a path to a log file or 
    # a data.frame (for testthat tests)
    if (is.character(log_file)) {
        log_table <- loadConsoleLog(log_file)
    } else {
        log_table <- log_file
    }

    task <- unique(log_table[['task']])

    if (length(task)>1) {
        stop("more than one task in the same log")
    }

    if (task == "ParseDb-split") {
        # step          task     field                    value
        #    1 ParseDb-split      FILE Sample8_quality-pass.tsv <- this input file
        #    1 ParseDb-split     FIELD               productive
        #    1 ParseDb-split NUM_SPLIT                     None
        #    1 ParseDb-split   OUTPUT1 Sample8_productive-F.tsv <- has two output files
        #    1 ParseDb-split   OUTPUT2 Sample8_productive-T.tsv <- 
        #    1 ParseDb-split   RECORDS                       92
        #    1 ParseDb-split     PARTS                        2
        #    1 ParseDb-split     SIZE1                        3
        #    1 ParseDb-split     SIZE2                       89
        file_field <- grepl("FILE",log_table[['field']])
        if (sum(file_field)>1) {
            stop("Can't deal with multiple input files.")
        }
        log_table[['field']][file_field] <- 'FILE1'
        records_field <- grepl("RECORDS",log_table[['field']])
        log_table[['field']][records_field] <- 'RECORDS1'
        # step          task     field                    value
        # 3     1 ParseDb-split     FILE1 Sample8_quality-pass.tsv <- file1
        # 4     1 ParseDb-split     FIELD               productive
        # 5     1 ParseDb-split NUM_SPLIT                     None
        # 6     1 ParseDb-split   OUTPUT1 Sample8_productive-F.tsv
        # 7     1 ParseDb-split   OUTPUT2 Sample8_productive-T.tsv
        # 8     1 ParseDb-split  RECORDS1                       92 <- has size records1
        # 9     1 ParseDb-split     PARTS                        2
        # 10    1 ParseDb-split     SIZE1                        3
        # 11    1 ParseDb-split     SIZE2                       89
        
        parts <- as.numeric(log_table[['value']][log_table[['field']]=="PARTS"])
        i <- 2
        while (i <= parts) {
            log_table <- bind_rows(log_table,
                      log_table[which(file_field),]
            )
            log_table[['field']][nrow(log_table)] <- paste0("FILE",i)
            log_table <- bind_rows(log_table,
                                   log_table[which(records_field),]
            )     
            log_table[['field']][nrow(log_table)] <- paste0("RECORDS",i)
            i <- i + 1
        }
        # df is now formatted for future pivot. inputs and outputs are
        # paired by number (eg file1 output1) and similar for
        # input records processed and size of the output files.
        # step          task     field                    value
        #    1 ParseDb-split     FILE1 Sample8_quality-pass.tsv
        #    1 ParseDb-split     FIELD               productive
        #    1 ParseDb-split NUM_SPLIT                     None
        #    1 ParseDb-split   OUTPUT1 Sample8_productive-F.tsv
        #    1 ParseDb-split   OUTPUT2 Sample8_productive-T.tsv
        #    1 ParseDb-split  RECORDS1                       92
        #    1 ParseDb-split     PARTS                        2
        #    1 ParseDb-split     SIZE1                        3
        #    1 ParseDb-split     SIZE2                       89
        #    1 ParseDb-split     FILE2 Sample8_quality-pass.tsv <- copied for output2
        #    1 ParseDb-split  RECORDS2                       92 <- copied for output2
    }
  
    if (task == "MakeDB-igblast") {
        log_table <- log_table %>%
            filter(field != "SEQ_FILE")
    }
    # Add RECORDS field to track input file size
    # TODO: update for multi input/output
    if ( "RECORDS" %in% log_table[["field"]] == FALSE ) {
        records <- log_table %>%
            tidyr::extract(field, 
                           into = c("field", "input_id"), 
                           regex = "([A-Za-z]+)([0-9]*)$", remove = FALSE) %>%
            group_by(input_id) %>%
            filter(field %in% c("FAIL","PASS")) %>%
            summarize(value=as.character(sum(as.numeric(value)))) %>%
            mutate(field = paste0("RECORDS", input_id),
                   step=log_table[["step"]][1],
                   task=log_table[["task"]][1])
        
        log_table <- bind_rows(
            log_table,
            records %>% select(-input_id)
        )
        
    }
    
    .getFileID <- function(x) {
        sapply(x, function(i){
            if (grepl("[0-9]$",i)){
                gsub("(.*?)([0-9]+)$","\\2",i)
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
                                values_from=value, names_from=field_type),
                by="file_id"
            )
    }
    log_table_input <- log_table_input %>%
        rename(input=value)

    log_table_output <- log_table %>%
        filter(!grepl("input", field_type)) %>%
        pivot_wider(id_cols=file_id,
                    values_from=value, names_from=field_type)

    log_table <- log_table_input %>%
        left_join(log_table_output, by="file_id")

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

#' @export
consoleLogsAsGraphs <- function(logs, metadata=NULL) {
    
    if (!is.null(metadata)) {
        if (nrow(metadata)>1) {
            # Address duplicated files coming from different originating samples
            # If needed, make filenames unique by adding sample_id. 
            # Relevant for 10x studies, where different samples can have files 
            # with the same filename (airr_rearrangement.tsv).
            .makeFile <- function(.data) { paste(.data[["sample_id"]], .data[["filename"]], sep=": ")}
            metadata <- metadata %>%
                group_by(filename) %>%
                mutate(n=n()) %>%
                ungroup() %>%
                rowwise() %>%
                mutate(name = ifelse(n>1, .makeFile(.data) , filename)) %>%
                select(-n)
        } else {
            stop("`metadata` is empty. Please update it or use NULL.")
        }
    }
    
    # .makeVertexName is a helper function to manage duplicated inputs in the logs. 
    # This is to allow for duplicated input files.
    # We found situations where the basename is duplicated i.e. multiple folders,
    # one per sample, and each folder has an airr_rearrangement.tsv file.
    # The log from the rename step in nf-core/airrflow only records the file name
    # (basename), because that is how Nextflow's path() works.
    # .makeUnique will add the sample id to the "duplicated" input files in 
    # the RenameFile tasks. graph_from_data_frame won't work with duplicated node
    # names.
    # .makeVertexName <- function(.data) {
    #     if (.data[["task"]] == "RenameFile") {
    #         sample_id <- sub("\\.[^\\.]*$","",.data[["output"]])
    #         paste0(sample_id,": ", .data[["input_id"]], collapse="")
    #     } else {
    #         # We should expect duplicated names only in the initial RenameFile tasks
    #         # as downstream processes add modifiers that will make file names unique (sample_id)
    #         stop("Unexpected duplicated input names.")
    #     }
    # }
    
    # Concatenate input and output with their sizes to create the graph nodes.
    # This should help to avoid issues with duplicated file names
    logs <-  logs %>% 
        rowwise() %>%
        mutate(input_id=paste(input,input_size, sep="_"),
               output_id=paste(output,output_size, sep="_")
               )
    # %>%
    #     group_by(input_id) %>%
    #     mutate(group_id=cur_group_id(),
    #            n=n(),
    #            is_dup=n>1 & length(unique(log_id))>1,
    #            n_idx=1:n()) %>%
    #     ungroup()
# 
#     logs <- logs %>%
#         rowwise() %>%
#         mutate(name = ifelse(is_dup, .makeVertexName(.data), input_id))
    
    # Create nodes from the input_id and output_id (uses file name a file size)
    # Remove duplicated nodes. This happens because outputs of some tasks are
    # inputs for others tasks, and they have the same size.
    vertex_meta <- bind_rows(
        logs %>% 
            select(input_id, input, input_size) %>% 
            rename(name=input_id) %>%
            rename(num_seqs=input_size, filename=input),
        logs %>% 
            select(output_id, output, output_size) %>%
            rename(name=output_id) %>%
            rename(num_seqs=output_size,filename=output),
    ) %>%
    distinct() 
    
    dup_v <- duplicated(vertex_meta$name)
    if (any(dup_v)) {
        stop("Unexpected duplicated vertex names.")
    }
    g <- graph_from_data_frame(logs %>% 
                                   select(input_id, output_id, task) %>%
                                   rename(name=input_id),
                               directed = TRUE,
                               vertices = vertex_meta)
    
    # Identify graph components, belonging to the different
    # groups of files that belong to the same original data
    # source
    V(g)$component <- igraph::components(g)$membership
    
    # Label the root nodes, the original input files
    g <- igraph::set_vertex_attr(g, "is_input", index = V(g), FALSE)
    g <- igraph::set_vertex_attr(g, "is_input", index = V(g)[degree(g, mode="in")==0], TRUE)
    
    # Add source task as node attribute
    tasks <- as_long_data_frame(g) %>%
        select(to_name, task)
    g <- igraph::set_vertex_attr(g, "task", index = V(g), "input")
    for (i in 1:nrow(tasks)) {
        this_task <- tasks[["task"]][i]
        this_vertex_name <- tasks[["to_name"]][i]
        v_idx <- V(g)$name == this_vertex_name
        g <- igraph::set_vertex_attr(g, "task", index = V(g)[v_idx], this_task)
    }
    
    # Propagate metadata and root name through the component  
    if (!is.null(metadata)) {
        for (i in unique(V(g)$component)) {
            c_idx <- V(g)$component == i
            this_meta <- metadata[metadata$name %in% V(g)$name[c_idx],,drop=FALSE] %>%
                select(-filename, -name)
            if (nrow(this_meta)==0) {
                fasta_idx <- grepl("\\.fasta_[0-9]+$|\\.fasta$", V(g)$name[c_idx])
                if (any(fasta_idx)) {
                    v_sample_id <- sub("([^_]+)_.+","\\1", V(g)$name[c_idx][fasta_idx])
                    this_meta <- metadata[metadata$sample_id %in% v_sample_id,,drop=FALSE] %>%
                        select(-filename, -name) %>%
                        distinct()
                }
            }
            this_root_name <- igraph::vertex_attr(g, "filename")[igraph::vertex_attr(g, "is_input") & c_idx]
            g <- igraph::set_vertex_attr(g, "root_name", index = V(g)[c_idx], this_root_name)
            if (nrow(this_meta)!=1 ) { stop("Expecting one row with metadata for this component.") }
            for ( meta_field in colnames(this_meta)) {
                g <- igraph::set_vertex_attr(g, meta_field, index = V(g)[c_idx], this_meta[[meta_field]])
            }
        }
    }
    
    # Example structure. Note rows 12 and 15.
    # name is the vertex name, that must be unique.
    # filename is the original filename
    # is_input labels the original input files
    #
    # > V(g)[[]]
    # + 20/20 vertices, named, from b2de37a:
    #                                name                   filename num_seqs component is_input                task                        root_name sample_id
    # 1         sample_x_productive-T.tsv  sample_x_productive-T.tsv       79         1    FALSE       ParseDb-split sample_x: airr_rearrangement.tsv  sample_x
    # 2         sample_y_quality-pass.tsv  sample_y_quality-pass.tsv       87         2    FALSE       FilterQuality sample_y: airr_rearrangement.tsv  sample_y
    # 3                      sample_x.tsv               sample_x.tsv       79         1    FALSE          RenameFile sample_x: airr_rearrangement.tsv  sample_x
    # 4        sample_y_junction-pass.tsv sample_y_junction-pass.tsv       87         2    FALSE  FilterJunctionMod3 sample_y: airr_rearrangement.tsv  sample_y
    # 5         sample_x_quality-pass.tsv  sample_x_quality-pass.tsv       79         1    FALSE       FilterQuality sample_x: airr_rearrangement.tsv  sample_x
    # 6              sample_y_db-pass.tsv       sample_y_db-pass.tsv       87         2    FALSE      MakeDB-igblast sample_y: airr_rearrangement.tsv  sample_y
    # 7          sample_y_sequences.fasta   sample_y_sequences.fasta       87         2    FALSE     ConvertDb-fasta sample_y: airr_rearrangement.tsv  sample_y
    # 8            sample_x_meta-pass.tsv     sample_x_meta-pass.tsv       79         1    FALSE         AddMetadata sample_x: airr_rearrangement.tsv  sample_x
    # 9            sample_y_meta-pass.tsv     sample_y_meta-pass.tsv       87         2    FALSE         AddMetadata sample_y: airr_rearrangement.tsv  sample_y
    # 10       sample_x_junction-pass.tsv sample_x_junction-pass.tsv       79         1    FALSE  FilterJunctionMod3 sample_x: airr_rearrangement.tsv  sample_x
    # 11            sample_x_igblast.fmt7      sample_x_igblast.fmt7       79         1    FALSE AssignGenes-igblast sample_x: airr_rearrangement.tsv  sample_x
    # 12 sample_y: airr_rearrangement.tsv     airr_rearrangement.tsv       87         2     TRUE               input sample_y: airr_rearrangement.tsv  sample_y
    # 13                     sample_y.tsv               sample_y.tsv       87         2    FALSE          RenameFile sample_y: airr_rearrangement.tsv  sample_y
    # 14        sample_y_productive-T.tsv  sample_y_productive-T.tsv       87         2    FALSE       ParseDb-split sample_y: airr_rearrangement.tsv  sample_y
    # 15 sample_x: airr_rearrangement.tsv     airr_rearrangement.tsv       79         1     TRUE               input sample_x: airr_rearrangement.tsv  sample_x
    # 16             sample_x_db-pass.tsv       sample_x_db-pass.tsv       79         1    FALSE      MakeDB-igblast sample_x: airr_rearrangement.tsv  sample_x
    # 17         sample_x_sequences.fasta   sample_x_sequences.fasta       79         1    FALSE     ConvertDb-fasta sample_x: airr_rearrangement.tsv  sample_x
    # 18            sample_y_igblast.fmt7      sample_y_igblast.fmt7       87         2    FALSE AssignGenes-igblast sample_y: airr_rearrangement.tsv  sample_y
    # 19          sample_x__scqc-pass.tsv    sample_x__scqc-pass.tsv       78         1    FALSE        SingleCellQC sample_x: airr_rearrangement.tsv  sample_x
    # 20          sample_y__scqc-pass.tsv    sample_y__scqc-pass.tsv       82         2    FALSE        SingleCellQC sample_y: airr_rearrangement.tsv  sample_y
    
    g_decomposed <- decompose(g)
    # return named list
    names(g_decomposed) <- sapply(g_decomposed, function(g) {
        if ("sample_id" %in% igraph::vertex_attr_names(g)) {
            igraph::vertex_attr(g, "sample_id")[igraph::vertex_attr(g, "is_input")]
        } else {
            igraph::vertex_attr(g, "name")[igraph::vertex_attr(g, "is_input")]
        }
    })
    list(
        "by_sample"=g_decomposed,
        "workflow"=g
    )
}


#' @export
cleanPlotData <- function(pd) {
    pd %>%
        group_by(component) %>%
        mutate(step = abs(y-max(y))) %>%
        arrange(step,x) %>%
        ungroup() %>%
        select(!any_of(c("x", "y", "component", "circular"))) %>%
        select(!starts_with(".")) %>%
        select(step, task , name, filename, num_seqs, everything()) %>%
        mutate(name = factor(name, levels=unique(name), ordered = T))
}

#' @export
plotLogGraph <- function(g) {
    gp <- ggraph(g, layout="tree") +
        #geom_node_point(aes(size = num_seqs), colour = "black") +
        geom_node_label(aes(label=paste0(name,": ", num_seqs, " sequences")),
                        colour = "black", label.size=unit(0.25, 'mm')) +
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
    gp
}

#' @export
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
    p <- ggraph(g3, layout = "tree") +
    # geom_edge_link0(aes(col=1,edge_width=2), arrow = arrow(length = unit(8, 'mm')))+
        geom_edge_link(color="darkblue", edge_width=unit(1, 'mm'), 
                       arrow = arrow(type = "closed", length = unit(3, 'mm')), 
                       end_cap = circle(3, 'mm'))+
        geom_node_point(shape=21,col="white",fill="black",size=5,stroke=0.5) +
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
#' @param log_files Vector of paths to log files
#' @export
readConsoleLogs <- function(log_files) {
    # log_files <-  strsplit(log_files,"[ ,]")[[1]]
    logs_list <- lapply(log_files, enchantr:::formatConsoleLog)
    names(logs_list) <- paste0("log_", 1:length(logs_list))
    bind_rows(logs_list, .id="log_id")
}


#' @param log_df data.frame of log data
#' @param style  Style plot. \code{workflow} will plot one figure with
#'               the pipeline steps. \code{decompose} will plot one figure
#'               for each starting input file.
#' @seealso  See also \link{readConsoleLogs}.
#' @export
plotConsoleLogs <- function(log_df, style=c("decompose", "workflow"), metadata=NULL) {
    style <- match.arg(style)
    logs <- consoleLogsAsGraphs(log_df, metadata)
    if (style == "decompose") {
        logs$by_sample
    } else {
        logs$workflow
    }

}

#' @export
graphDataFrame <- function(g) {
    exclude_cols <- c("from", "to",
                      "from_name", "from_task", "from_component",
                      "to_task")
    .fix_from_names <- function(x) {
        sub("^from_","",x)
    }
    as_long_data_frame(g) %>%
        rename(input=from_filename,
               output=to_filename,
               task=task,
               input_size=from_num_seqs,
               output_size=to_num_seqs,
               is_input=from_is_input) %>%
        select(!any_of(exclude_cols)) %>%
        select(!starts_with("to_")) %>%
        rename_with(.fix_from_names) %>%
        select(any_of(c("sample_id", "input", "input_size", 
                        "task", "output", "output_size")),
               everything())
}
