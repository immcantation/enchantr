#' Create an Immcantation ...
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
collapse_duplicates_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "collapse_duplicates_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 


#' @export
findDuplicates <- function (db, groups="sample_id",
                            id = "sequence_id", 
                            seq = "sequence_alignment", 
                            text_fields = NULL, 
                            num_fields = c("conscount", "dupcount"), 
                            seq_fields = NULL,
                            add_count = TRUE,
                            ignore = c("N", "-", ".", "?"), sep=",",
                            dry = F, verbose = F,
                            nproc=1) {
    
    c_call <- NULL
    if ("c_call" %in% colnames(db)) { 
        if (any(!is.na(db[['c_call']]))) {
            c_call <- "c_call"   
        }
    }
    columns <- c(groups, id, seq, text_fields, num_fields, seq_fields, 
                 "v_call", "d_call", "j_call", "junction_length", c_call, "productive")
    columns <- columns[!is.null(columns)]
    check <- alakazam:::checkColumns(db, columns)
    if (!check == TRUE ) { stop(check) }
    
    db[['row_idx']] <- 1:nrow(db)
    
    db[['collapse_idx']] <- db %>%
        mutate(v_gene=getGene(v_call),
               d_gene=getGene(d_call),
               j_gene=getGene(j_call)) %>%
        group_by(!!!rlang::syms(c(groups, "v_gene", "j_gene", c_call, "junction_length", "productive"))) %>%
        group_indices()
    
    
    # If the user has previously set the cluster and does not wish to reset it
    if(!is.numeric(nproc)){
        stop_cluster <- FALSE
        cluster <- nproc
        nproc <- 0
    } else {
        stop_cluster <- TRUE
    }
    # Ensure that the nproc does not exceed the number of cores/CPUs available
    nproc <- min(nproc, cpuCount())
    
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters &
    # export all nesseary environment variables, functions and packages.
    if (nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else {
        cluster_type <- "FORK"
        if (.Platform$OS.type == "windows") {
            cluster_type <- "PSOCK"
        }
        if (nproc != 0) {
            #cluster <- makeCluster(nproc, type="SOCK")
            cluster <- parallel::makeCluster(nproc, type= cluster_type)
        }
        if (cluster_type == "PSOCK") {
            parallel::clusterExport(cluster,
                                    list('db', id, seq, text_fields,
                                         num_fields, seq_fields, add_count, 
                                         ignore, sep, dry, verbose, columns),
                                    envir=environment() )
        }
        registerDoParallel(cluster)
    }
    
    group_idx <- unique(db[['collapse_idx']])
    
    collapse_pass <- bind_rows(foreach(i=1:length(group_idx),
                             .verbose=FALSE,
                             .errorhandling='stop') %dopar% {
                              
                              this_group <- group_idx[i]
                              
                              collapsed_db <- collapseDuplicates(db %>%
                                                                     filter(collapse_idx == this_group),
                                                     id = id,
                                                     seq = seq,
                                                     text_fields = text_fields, 
                                                     num_fields = num_fields, 
                                                     seq_fields = seq_fields,
                                                     add_count = add_count,
                                                     ignore = ignore, sep=sep,
                                                     dry = dry, verbose = verbose)
                              if (add_count) {
                                  columns <- c(columns, 'collapse_count','row_idx')
                              }
                              collapsed_db %>%
                                  select(!!!rlang::syms(columns))
                             })
    
    if (stop_cluster & !is.numeric(nproc)) {
        parallel::stopCluster(cluster)
    }
    
    db %>%
        mutate(collapse_pass=row_idx %in% collapse_pass[['row_idx']]) %>%
        arrange(row_idx) %>%
        select(-collapse_idx, -row_idx)
}
