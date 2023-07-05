#' Create an Immcantation ...
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
find_threshold_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "find_threshold_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 


# TODO: document, return
# TODO: example
# TODO: test
#' @seealso  See also \link{shazam::findThreshold}. 
#' @export
findThresholdDb <- function(db, distanceColumn="dist_nearest", 
                            crossDistanceColumn="cross_dist_nearest",
                            method=c("gmm", "dens"), 
                            model=c("gamma-gamma"),
                            edge = 0.9, subsample=NULL, 
                            cutoff = "user", spc = 0.995,
                            nproc=1, fields=NULL ) {
    # Hack for visibility of data.table and foreach index variables
    idx <- yidx <- .I <- NULL
    
    # Initial checks
    method <- match.arg(method)
    
    if (!is.data.frame(db)) { stop('Must submit a data frame') }
    # Check for valid columns
    columns <- c(distanceColumn, fields)
    columns <- columns[!is.null(columns)]
    
    check <- alakazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }    
    
    if (!is.null(crossDistanceColumn)) {
        if (!crossDistanceColumn %in% colnames(db)) {
            stop("The column ",crossDistanceColumn," does not exist in `db`")
        } 
    }
    # Convert the db (data.frame) to a data.table & set keys
    # This is an efficient way to get the groups of V J L, instead of doing dplyr
    dt <- data.table(db)
    
    # Get the group indexes
    dt <- dt[, list( yidx = list(.I) ) , by = fields ]
    groups <- dt[,yidx]
    lenGroups <- length(groups) 
    
    # Create cluster of nproc size and export namespaces
    # If user wants to paralellize this function and specifies nproc > 1, then
    # initialize and register slave R processes/clusters & 
    # export all nesseary environment variables, functions and packages.
    if( nproc==1 ) {
        # If needed to run on a single core/cpu then, register DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # Export groups to the clusters
    if (nproc > 1) { 
        export_functions <- list("db",
                                 "groups", 
                                 "distanceColumn", "crossDistanceColumn",
                                 "findThreshold",
                                 "method", "edge", "subsample", "cutoff", "spc",
                                 "plotGmmThreshold")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    } 
    
    list_db <- foreach(idx=iterators::icount(lenGroups), .errorhandling='stop') %dopar% {
        db_group <- db[groups[[idx]], ]
        if (any(!is.na(db_group[[distanceColumn]]))) {
            if (!is.null(crossDistanceColumn)) {
                cross_distances <- db_group[[crossDistanceColumn]]
                if (all(is.na(cross_distances))) {
                    warning("All 'cross_dist_nearest` values are NA, will use cross='NULL'.")
                    cross_distances <- NULL
                }
            } else {
                cross_distances <- NULL
            }
            threshold <- findThreshold(db_group[[distanceColumn]], 
                                       cross=cross_distances,
                                       method=method, 
                                       model=model,
                                       cutoff = cutoff, spc = spc,
                                       edge=edge, subsample=subsample)
        
            # threshold_vector <- sapply(slotNames(threshold), function(this_slot) {
            #     slot(threshold, this_slot)
            # })
            # threshold_vector <- threshold_vector[names(threshold_vector)!="x"]
            # if (!is.null(fields)) {
            #     c(unique(db_group[,fields]),threshold_vector)                
            # } else {
            #     threshold_vector
            # }
            p <- NULL
            if (!is.null(threshold)) {
                p <- plotGmmThreshold(threshold, cross=cross_distances, silent=TRUE)
            }
        } else {
            threshold <- NULL
            p <- NULL
        }
        list("fields"=unique(db_group[,fields]),
             "GmmThreshold"=unclass(threshold),
             "plot"=p
        )
        
    } 
    
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    #df <- dplyr::bind_rows(list_db)
    #df
    
    list_db
}

#' @export
gmmSummary <- function(gmm) {
    l <- lapply(gmm, function(thisGmm) {
        fields <- as.data.frame(thisGmm[['fields']], nrows=1)
        colnames(fields) <- paste0("fields", c("",1:(ncol(fields))))[1:ncol(fields)]
        gmmFit <- thisGmm[['GmmThreshold']]
        if (!is.logical(gmmFit)) {
            slots <- slotNames(getClass("GmmThreshold"))[slotNames(getClass("GmmThreshold")) != "x"]
            if (!is.null(gmmFit)) {
                fit <- sapply(slots, function(s) {
                    slot(gmmFit, s)
                })
            } else {
                fit <- rep(NA, length(slots))
                names(fit) <- slots
            }
            
        } else {
            fit <- data.frame()
        }
        c(fields, fit)
    })
    df <- bind_rows(l)
    # modify df for easier summary viewing
    df$threshold <- as.numeric(unlist(df$threshold))
    df <- subset(df, select = -c(a1, a2, b1, b2, c1, c2) )
    df
}
