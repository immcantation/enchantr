#' Create an Immcantation ...
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
collapse_duplicates_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                                "templates", "project",
                                "collapse_duplicates_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}


#' @export
findDuplicates <- function (db, 
                            groups="sample_id",
                            id = "sequence_id",
                            seq = "sequence_alignment",
                            text_fields = NULL,
                            num_fields = c("consensus_count", "duplicate_count"),
                            seq_fields = NULL,
                            add_count = TRUE,
                            ignore = c("N", "-", ".", "?"), 
                            sep=",",
                            dry = F, verbose = F,
                            nproc=1) {
    
    # Including c_call if present in the groups, so sequences with different c_call will not be collapsed
    c_call <- NULL
    if ("c_call" %in% colnames(db)) {
        if (any(!is.na(db[['c_call']]))) {
            c_call <- "c_call"
        }
    }

    d_call <- NULL
    if ("d_call" %in% colnames(db)) {
        if (any(!is.na(db[['d_call']]))) {
            d_call <- "d_call"
        }
    }

    # For many bulk protocols, the isotype will be determined with either the cprimer or the cregion annotations
    # Including cprimer and cregion if present in the groups, so sequences with different cprimer or cregions will not be collapsed
    c_primer <- NULL
    if ("c_primer" %in% colnames(db)) {
        if (any(!is.na(db[['c_primer']]))) {
            c_primer <- "c_primer"
        }
    }

    cregion <- NULL
    if ("cregion" %in% colnames(db)) {
        if (any(!is.na(db[['cregion']]))) {
            cregion <- "cregion"
        }
    }

    columns <- c(groups, id, seq, text_fields, num_fields, seq_fields,
                "v_call", d_call, "j_call", "junction_length", c_call, "productive", c_primer, cregion)
    columns <- columns[!is.null(columns)]
    check <- alakazam::checkColumns(db, columns)
    if (!check == TRUE ) { stop(check) }

    db[['finddups_row_idx']] <- 1:nrow(db)
    db[['seq_len']] <- sapply(db[[seq]], nchar)

    db[['collapse_idx']] <- db %>%
        mutate(v_gene=getGene(v_call),
                d_gene=getGene(d_call),
                j_gene=getGene(j_call)) %>%
        group_by(!!!rlang::syms(c(groups, "v_gene", "j_gene", c_call, c_primer, cregion, "junction_length", "productive", "seq_len"))) %>%
        group_indices()

    db_subset <- db %>%
        select(all_of(c(columns, "finddups_row_idx", "collapse_idx")))
    db <- db %>%
        select(!any_of(c(columns, "collapse_idx")))

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
            logfile <- tempfile("log", fileext = ".txt")
            message("Cluster log file: ", logfile)
            cluster <- parallel::makeCluster(nproc, type= cluster_type, outfile = logfile)
        }
        if (cluster_type == "PSOCK") {
            parallel::clusterExport(cluster,
                                    list('db_subset', 'id', 'seq', 'text_fields',
                                            'num_fields', 'seq_fields', 'add_count',
                                            'ignore', 'sep', 'dry', 'verbose', 'columns'),
                                    envir=environment() )
        }
        registerDoParallel(cluster)
    }

    group_idx <- unique(db_subset[['collapse_idx']])

    collapse_pass <- bind_rows(
        foreach(i=iterators::icount(length(group_idx)), .verbose=FALSE, .errorhandling='stop') %dopar% {

            this_group <- group_idx[i]
            this_group_size <- sum(db_subset[["collapse_idx"]] == this_group)
            if (this_group_size == 1 ) {
                return(
                    db_subset %>%
                    filter(collapse_idx == this_group)
                    )
                }

            collapsed_db <- collapseDuplicates(db_subset %>%
                                                filter(collapse_idx == this_group),
                                                id = id,
                                                seq = seq,
                                                text_fields = text_fields,
                                                num_fields = num_fields,
                                                seq_fields = seq_fields,
                                                add_count = add_count,
                                                ignore = ignore,
                                                sep=sep,
                                                dry = dry, verbose = verbose)
            collapsed_db %>%
                select(any_of(c(columns,'collapse_count', 'finddups_row_idx')))
        },
    )

    if (stop_cluster & !is.numeric(nproc)) {
        parallel::stopCluster(cluster)
    }

    db %>%
        left_join(db_subset, by="finddups_row_idx") %>%
        mutate(collapse_pass=finddups_row_idx %in% collapse_pass[['finddups_row_idx']]) %>%
        left_join(collapse_pass %>% select(any_of(c("finddups_row_idx", "collapse_count"))), by="finddups_row_idx") %>%
        arrange(finddups_row_idx) %>%
        select(!any_of(c('collapse_idx', 'finddups_row_idx')))
}
