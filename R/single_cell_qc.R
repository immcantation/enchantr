#' Create an Immcantation ...
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
single_cell_qc_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package="enchantr"),"rstudio", "templates", "project", "single_cell_qc_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir,full.names = T)
    file.copy(project_files, project_dir, recursive=TRUE)
} 


# TODO: document
#' @export
countSequencesPerCell <- function(db) {
    seqs_per_cell <- db %>%
        group_by(sample_id, cell_id, locus) %>%
        summarize(cell_num_sequences=n(), 
                  cell_num_isotypes=length(unique(c_call)),
                  cell_isotypes=paste(unique(c_call),sep=",",collapse=",")) %>%
        ungroup() %>%
        arrange(sample_id, desc(cell_num_sequences))
}

# TODO: document
#' @export
plotSequencesPerCell <- function(seqs_per_cell) {
    ggplot(seqs_per_cell,aes(x=cell_num_sequences)) +
        geom_histogram(binwidth = 1) +
        scale_x_continuous(breaks=1:max(seqs_per_cell$cell_num_sequences)) +
        facet_wrap(sample_id~locus) +
        labs(
            x= "Number of sequences in a cell.",
            y=" Count",
            caption="Distribution of the number of sequences per cell and locus."
        )
}

#' @export
scQC <- function(db, cell_id="cell_id") {
    db[['chain']] <- getChain(db[["locus"]])
    db[['tmp_scqc_row']] <- 1:nrow(db)
    db <- db %>%
        group_by(sample_id, cell_id) %>%
        mutate(paired_chain= all(c("VL", "VH") %in% chain)) %>%
        ungroup() %>%
        group_by(sample_id, chain, cell_id) %>%
        mutate(num_sequences=n(),
               num_isotypes=length(unique(c_call)),
               cell_isotypes=paste(unique(c_call),sep=",",collapse=","),
               max_count=max(consensus_count),
               perc_of_cell_cons_count=100*consensus_count/sum(consensus_count),
               is_most_abundant= consensus_count==max(consensus_count) & sum(consensus_count==max(consensus_count)) == 1,
               cons_count_ratio=consensus_count/max_count) %>%
        group_by(sample_id, chain, cell_id) %>%   
        do(cellQC(.)) %>%
        ungroup() 
    qclog <- db %>%
        arrange(tmp_scqc_row) %>%
        select(
               sample_id, 
               cell_id, 
               sequence_id, 
               paired_chain,
               num_sequences, 
               num_isotypes, 
               cell_isotypes, 
               max_count, 
               perc_of_cell_cons_count, 
               is_most_abundant, 
               cons_count_ratio, 
               scqc_pass)
    list(
        "log"= qclog,
        "pass"=sum(qclog$scqc_pass),
        "fail"=sum(!qclog$scqc_pass)
    )
}

cellQC <- function(df) {
    if (length(unique(df$sample_id)) > 1) { stop ("Expectig data from one sample_id. Group by sample_id.")}
    if (length(unique(df$chain)) > 1) { stop ("Expectig data from one chain. Group by chain.")}
    df$scqc_pass <- FALSE
    if (nrow(df)==1) {
        df$scqc_pass <- T
    } else {
        df <- df %>%
            arrange(desc(perc_of_cell_cons_count))
        if (df$perc_of_cell_cons_count[1] >50) {
            df$scqc_pass[1] <- TRUE
        } 
        if (all(grepl("[LKlk]",df$chain)))  {
            if (sum(df$perc_of_cell_cons_count[1:2]) > 75 & df$perc_of_cell_cons_count[2] > 25 ) {
                df$scqc_pass[1:2] <- TRUE
            }
        }
    }
    if (sum(df$scqc_pass) == 0) {
        warning("No sequences selected for cell ", unique(df$cell_ide), "consensus_count: ", unique(df$consensus_count))
    }
    df %>%
        arrange(tmp_scqc_row)
}


#' findSingleCellDuplicates(db, "sample_id")
#' @export
findSingleCellDuplicates <- function(db, groups, 
                                     cell = "cell_id", 
                                     seq="sequence_alignment",
                                     sequence_id="sequence_id") {
    # Check for valid columns
    if (is.null(groups)) { stop("`groups` must be a valid column name")}
    columns <- c(groups,cell, seq, sequence_id)
    columns <- columns[!is.null(columns)]
    check <- alakazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Check that sequence_id are unique, because I will use this id
    # later to remove the duplicated sequences
    if (any(duplicated(db[[sequence_id]]))) {
        stop("Duplicated `sequence_id` found. Expecting unique sequence identifiers.")
    }
    
    # Identify groups
    db[['group_idx']] <- db %>%
        dplyr::group_by(!!rlang::sym(groups)) %>% 
        group_indices()
    
    # Identify cells present in more than one group
    db_shared <- db %>%
        select(!!!rlang::syms(c(columns, "group_idx"))) %>%
        group_by(!!rlang::sym(cell)) %>%
        mutate(is_shared_cell = length(unique(group_idx))>1) %>%
        ungroup() %>%
        filter(is_shared_cell)
    
    # Helper function to label as duplicate T/F those sequences
    # that have a 0 distance sequence in other group
    .isDuplicate <- function(x, seq) {
        sequences <- x[[seq]]
        
        # is_duplicate defaults to FALSE
        is_duplicate <- rep(F,length(sequences))
        
        # Calculate distance between sequences
        # use gap=0 to treat gaps as Ns
        # and label as T those with a distance value of 0
        pwd <- pairwiseDist(sequences, dist_mat=getDNAMatrix(gap=0)) == 0
        # don't use self comparison
        diag(pwd) <- NA
        # Get the indices for the duplicated sequences
        duplicated_idx <- which(pwd, arr.ind = T)
        if (nrow(duplicated_idx)>0) {
            # Group 
            dup_rows_group <- x[['group_idx']][duplicated_idx[,"row"]]
            dup_cols_group <- x[['group_idx']][duplicated_idx[,"col"]]
            dup_sequences <- duplicated_idx[,"row"][dup_rows_group != dup_cols_group]
            if(length(dup_sequences)>0) {
                is_duplicate[dup_sequences] <- T
            }
            
        }
        x$is_duplicate <- is_duplicate
        x
    }
    
    duplicated_seq_id <- db_shared %>%
        group_by(!!rlang::sym(cell)) %>%
        do(.isDuplicate(., seq)) %>%
        filter(is_duplicate) %>%
        select(-group_idx) %>%
        pull(!!rlang::sym(sequence_id))
    
    db[['sc_duplicate']] <- F
    if (length(duplicated_seq_id)>0) {
        db[['sc_duplicate']][db[[sequence_id]] %in% duplicated_seq_id] <- T
    }
    db %>%
        select(-group_idx)
}

#' removeSingleCellDuplicates(db, "sample_id")
#' @export
removeSingleCellDuplicates <- function(db, groups, 
                                       cell = "cell_id", 
                                       seq="sequence_alignment",
                                       sequence_id="sequence_id") {
    findSingleCellDuplicates(db, groups, cell, seq, sequence_id) %>%
        filter(sc_duplicate == F) %>%
        select(-sc_duplicate)
}



##---------------------- TMP -------------------------- ##
# This function modified previous overlap heatmap function and is used to show, between two samples 
# in single cell data, the overlap of cells which has same barcode and sequences
# find cells having same seq and cell barcode
#' @export
singlecell_sharing_matrix <- function(db, 
                                      groups='sample_id', 
                                      cell = "cell_id",
                                      seq="sequence_alignment",
                                      order_by=NULL
) {
    
    # Check for valid columns
    if (is.null(groups)) { stop("`groups` must be a valid column name")}
    columns <- c(groups,cell, seq, order_by)
    columns <- columns[!is.null(columns)]
    check <- alakazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # Identify groups
    db[['group_idx']] <- db %>%
        dplyr::group_by(!!rlang::sym(groups)) %>% 
        group_indices()
    
    # extract group names
    groups_table = unique(db[,c("group_idx",groups), drop=F])
    
    groups_table$group_name <- sapply(1:nrow(groups_table), function(g) {
        use <- colnames(groups_table) %in% c(groups) == T
        paste(as.matrix(groups_table[g,use,drop=F]),collapse="_")
    } )
    
    num_groups = nrow(groups_table)
    
    # define clone overlap matrices to be plotted
    overlap_matrix <- matrix(0, num_groups, num_groups)        
    dimnames(overlap_matrix) <- list(groups_table$group_name, groups_table$group_name)
    
    overlap_count <- overlap_total <- overlap_percent <- overlap_matrix # matrix of overlap counts

    # calculate the number of overlapped cells
    for (ii in 1:nrow(groups_table)) {    
        for (jj in 1:nrow(groups_table)) {
            
            ii_group_idx <- groups_table[['group_idx']][ii]
            jj_group_idx <- groups_table[['group_idx']][jj]
            
            # cell in groups i and j
            cells_ii = db[[cell]][db[['group_idx']] == ii]
            cells_jj = db[[cell]][db[['group_idx']] == jj]
            # seqs in groups i and j
            seqs_ii = db[[seq]][db[['group_idx']] == ii]
            seqs_jj = db[[seq]][db[['group_idx']] == jj]
            
            # find the number of cells having same seq and cell barcode shared between groups i and j
            # TO reduce match time, find shared cells first and then compare seq for those cells
            shared_cells = unique(intersect(cells_ii, cells_jj))
            index_ii = match(shared_cells, cells_ii)
            index_jj = match(shared_cells, cells_jj)
            
            num_shared_cells = sum(seqs_ii[index_ii] == seqs_jj[index_jj])
            
            total_ii_cells = length(unique(cells_ii))
            total_jj_cells = length(unique(cells_jj))
            total_cells = min(total_ii_cells, total_jj_cells)
            #total_cells = length(union(seqs_ii,seqs_jj))
            
            # find the percentage of cells shared relative to the total number of cells
            perc_shared_cells = round (num_shared_cells / total_cells * 100 , 4)
            
            # Fill in matrices
            overlap_count[ii, jj] <- num_shared_cells
            overlap_total[ii, jj] <- total_cells
            overlap_percent[ii, jj] = perc_shared_cells
            
            #if(ii==jj){ overlap_perc_matrix[ii, jj] = length(seqs_ii) }
        }
    }
    if(!is.null(order_by)){
        overlap_count <- overlap_count[order_by,order_by]
        overlap_total <- overlap_total[order_by,order_by]
        overlap_percent <- overlap_percent[order_by,order_by]
    }
    return(createOverlap(overlap_count, overlap_total, overlap_percent))
    
}

# S4 class defining a Overlap object
setClass("Overlap", 
         slots = c( overlap_count="matrix",
                    overlap_total="matrix",
                    overlap_percent="matrix" )
)


createOverlap <- function( overlap_count=overlap_count,
                           overlap_total=overlap_total,
                           overlap_percent=overlap_percent ) {
    
    # Define RegionDefinition object
    overlap <- new("Overlap",
                   overlap_count=overlap_count,
                   overlap_total=overlap_total,
                   overlap_percent=overlap_percent)
    
    return(overlap)
}
# plotOverlapResubmission( list_objOverlap=list(uniqueSeq, uniqueClones), titleText="", 
#                          order_by=colnames(uniqueSeq@overlap_count), sizeCount=3 , subject=subject)
plotOverlapSingleCell <- function(list_objOverlap="", 
                                  titleText="", 
                                  ROI=NULL, 
                                  order_by=NULL,
                                  sizeDiag=3.5,
                                  sizeCount=5,
                                  sizePercent=4,
                                  A=NULL,
                                  B=NULL,
                                  subject=""){
    
    list_data_m <- list()
    
    for(x in 1:length(list_objOverlap)){
        
        objOverlap <- list_objOverlap[[x]]
        if(!is.null(ROI)){
            if(length(ROI)>1){
                roi <- match(ROI,row.names(objOverlap@overlap_count))
                roi <- roi[!is.na(roi)]
            }else{
                roi <- which(grepl(ROI,row.names(objOverlap@overlap_count)))
            }
            objOverlap@overlap_count <- objOverlap@overlap_count[roi,roi]
            objOverlap@overlap_total <- objOverlap@overlap_total[roi,roi]
            objOverlap@overlap_percent <- objOverlap@overlap_percent[roi,roi]
        }    
        
        
        data_m_count <- reshape::melt(objOverlap@overlap_count)
        data_m_total <- reshape::melt(objOverlap@overlap_total)
        data_m <- reshape::melt(objOverlap@overlap_percent)
        
        #data_m$Y1 <- cut(data_m$value,breaks = c(0,0.00001,seq(0.00002,99.99998,by=1),99.99999,100),right = FALSE,include.lowest=TRUE) 
        data_m$Y1 <- cut(data_m$value,breaks = c(0,0.00001,1,2,3,4,6,10,20,30,99.99999,100),right = FALSE,include.lowest=TRUE)  
        #copa <- colorRampPalette(brewer.pal(5,"Greys"))(length(levels(data_m$Y1)))
        copa <- colorRampPalette(brewer.pal(9,"Oranges"))(length(levels(data_m$Y1)))
        names(copa) <- levels(data_m$Y1)
        
        names(data_m)[3] <- "Percent"        
        data_m$Count <- data_m_count$value
        data_m <- data_m[ is.finite(data_m$Percent), ]
        
        #Format value to print
        data_m$Value = round(data_m$Percent,2)
        data_m$Value[ data_m$Value==0 ] = ""
        
        # Non diagonal are counts/percent
        #data_m$Value[ data_m$Value!="" ] = paste0( data_m$Count[ data_m$Value!="" ] , "\n(", data_m$Percent[ data_m$Value!="" ] , "%)")
        data_m$Value[ data_m$Value!="" ] = data_m$Count[ data_m$Value!="" ]
        data_m$Value[ as.vector(data_m$X1) == as.vector(data_m$X2) ] = ""
        
        data_m$ValuePercent <- data_m$Value
        data_m$ValuePercent[ data_m$Value!="" ] = data_m$Percent[ data_m$Value!="" ]
        
        data_m$ValuePercent <-
            unlist(
                sapply(data_m$ValuePercent, function(x){ 
                    if(x!="") {
                        return(paste0(format(round(as.numeric(x),1),nsmall=1),"%"))
                    } else {
                        return("")
                    }
                }))
        
        # Diagonals are counts
        data_m$DiagonalCount = ""
        data_m$DiagonalCount[as.vector(data_m$X1) == as.vector(data_m$X2)]= data_m$Count[as.vector(data_m$X1) == as.vector(data_m$X2)]
        
        data_m$X1 <- factor(data_m$X1, levels=order_by)
        data_m$X2 <- factor(data_m$X2, levels=order_by)
        
        list_data_m[[x]] <- data_m
    }
    
    
    data_m1 <- list_data_m[[1]]		
    data_m2 <- list_data_m[[2]]		
    
    data_m <- data_m1
    numbItems <- nrow(list_objOverlap[[1]]@overlap_count )
    pos <- as.vector(upper.tri(matrix(NA, nrow=numbItems,ncol=numbItems,byrow=T)))
    
    
    data_m$Value[pos] <-   data_m2$Value[pos]
    data_m$ValuePercent[pos] <-   data_m2$ValuePercent[pos]
    data_m$DiagonalCountLower <-   data_m$DiagonalCount
    data_m$DiagonalCount <-   data_m2$DiagonalCount
    data_m$Y1[pos] <-   data_m2$Y1[pos]
    #sizeDiag=3
    p <- ggplot(data=data_m, aes(X1, X2)) +
        geom_tile(aes(fill = Y1), colour = "white") +
        geom_text(aes(label = DiagonalCount), size=sizeDiag, angle = 50, color="white", vjust=-.5, hjust=0.5) +
        geom_text(aes(label = DiagonalCountLower), size=sizeDiag, angle = 50, color="white", vjust=1.5, hjust=0.5) +
        geom_text(aes(label = Value), size=sizeCount, color="black", vjust=0) +
        #geom_text(aes(label = ValuePercent), size=sizePercent, color="black", vjust=1.4) +
        scale_fill_manual(breaks=names(copa),values=copa, na.value="#FFFFFF", guide='none') +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black", size=10),
            axis.text.y = element_text(colour = "black",size=10,hjust = 0),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            text = element_text(size = 12, colour = "black"))+
        annotate("segment", x = 0.5, xend = numbItems+.5, y = 0.5, yend = numbItems+.5, colour = "white") + 
        ggtitle(titleText) + 
        labs(x="",y="")+
        theme(plot.title = element_text(lineheight=.8, face="bold"))+
        #theme(plot.margin = unit(c(0,0,.25,.25), "cm")) +
        theme(axis.title.y=element_text(vjust=1.3)) +
        theme(axis.title.x=element_text(vjust=-.4))
    
    return(p)
} 
