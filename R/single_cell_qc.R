#' Create an Immcantation Chimera detection
#' 
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation Input Validation
#' to create the skeleton of an Immcantation Input Validation project
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
        group_by(cell_id, locus) %>%
        summarize(cell_num_sequences=n(), 
                  cell_num_isotypes=length(unique(c_call)),
                  cell_isotypes=paste(unique(c_call),sep=",",collapse=",")) %>%
        arrange(desc(cell_num_sequences)) %>%
        ungroup() %>%
        arrange(desc(cell_num_sequences))
}

# TODO: document
#' @export
plotSequencesPerCell <- function(seqs_per_cell) {
    ggplot(seqs_per_cell,aes(x=cell_num_sequences)) +
        geom_histogram(binwidth = 1) +
        scale_x_continuous(breaks=1:max(seqs_per_cell$cell_num_sequences)) +
        facet_wrap(~locus) +
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
        group_by(cell_id) %>%
        mutate(paired_chain= all(c("VL", "VH") %in% chain)) %>%
        ungroup() %>%
        group_by(chain, cell_id) %>%
        mutate(num_sequences=n(),
               num_isotypes=length(unique(c_call)),
               cell_isotypes=paste(unique(c_call),sep=",",collapse=","),
               max_count=max(consensus_count),
               perc_of_cell_cons_count=100*consensus_count/sum(consensus_count),
               is_most_abundant= consensus_count==max(consensus_count) & sum(consensus_count==max(consensus_count)) == 1,
               cons_count_ratio=consensus_count/max_count) %>%
        group_by(chain, cell_id) %>%   
        do(cellQC(.)) %>%
        ungroup() %>%
        group_by(chain, cell_id) %>%
        arrange(desc(num_sequences), perc_of_cell_cons_count) %>%
        ungroup()
    qclog <- db %>%
        arrange(tmp_scqc_row) %>%
        select(cell_id, 
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