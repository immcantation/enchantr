#' Create an Immcantation ...
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
single_cell_qc_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                              "templates", "project",
                              "single_cell_qc_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }

    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}


#' countSequencesPerCell
#'
#' \code{countSequencesPerCell} counts the number of sequences of each isotype in each sample's cell
#'
#' @param    db           data.frame with AIRR-format style columns.
#' @param    sample_id    column in \code{db} containing sample identifiers
#' @param    cell_id      column in \code{db} containing cell identifiers
#' @param    locus        column in \code{db} containing locus assignments
#' @param    c_call       column in \code{db} containing constant region assignments
#'
#' @return   A data.frame with cell counts
#' @examples
#' data(Example10x, package="alakazam")
#' db <- Example10x
#' db[["sample_id"]] <- "example_sample"
#' countSequencesPerCell(db[1:10,])
#'
#' @seealso  See also \link{plotSequencesPerCell}.
#' @export
countSequencesPerCell <- function(db,
                                  sample_id="sample_id",
                                  cell_id="cell_id",
                                  locus="locus",
                                  c_call="c_call") {
    db %>%
      group_by(sample_id, cell_id, locus) %>%
      summarize(cell_num_sequences = n(),
                cell_num_isotypes = length(unique(c_call)),
                cell_isotypes = paste(unique(c_call), sep = ",", collapse = ","),
                .groups = "drop") %>%
      ungroup()
}

#' TODO: example
#' plotSequencesPerCell
#'
#' \code{plotSequencesPerCell} plots a histogram of the distribution of the
#' number of sequences per cell for each sample and locus.
#'
#' @param    seqs_per_cell  data.frame with the number of sequences per cell.
#'
#' @return   A ggplot
#'
#' @seealso  See also \link{countSequencesPerCell}.
#'
#' @export
plotSequencesPerCell <- function(seqs_per_cell) {
    ggplot(seqs_per_cell, aes(x = cell_num_sequences)) +
      geom_histogram(binwidth = 1, linewidth = 0.2, color = "black") +
      scale_x_continuous(breaks = 1:max(seqs_per_cell$cell_num_sequences)) +
      facet_grid(sample_id~locus) +
      labs(x = "Number of sequences per cell",
           y = "Number of cells",
           caption = "Distribution of the number of sequences per cell and locus")
}

#' TODO: example
#' findLightOnlyCells
#'
#' \code{findLightOnlyCells} identifies cells with only light chains.
#'
#' @param    db          data.frame with AIRR-format style columns.
#' @param    sample_id   column in \code{db} containing sample identifiers
#' @param    cell_id     column in \code{db} containing cell identifiers
#' @param    locus       column in \code{db} containing locus assignments
#' @param    fields      Columns in \code{db}, in addition to \code{sample_id},
#'                       that should be used to group sequences to be
#'                       analyzed independently.
#' @return   The input \code{db} with an additional column named \code{light_only_cell}
#'           with values TRUE/FALSE.
#' @export
findLightOnlyCells <- function(db,
                               sample_id = "sample_id",
                               cell_id="cell_id", locus="locus",
                               fields=NULL) {

    groups <- c(sample_id, fields)

    db <- db %>%
            group_by(!!!rlang::syms(c(cell_id, groups))) %>%
            mutate(light_only_cell = sum(isHeavyChain(!!rlang::sym(locus))) == 0) %>%
            ungroup()

    num_cells <- db %>%
                   dplyr::filter(light_only_cell) %>%
                   select(!!!rlang::syms(c(sample_id, cell_id, fields))) %>%
                   distinct() %>%
                   nrow()

    message("db size: ", nrow(db), " sequences.")
    message("number of light only cells: ", num_cells)

    db
}

#' TODO: example
#' removeDoublets
#'
#' \code{removeDoublets} removes cells with multiple heavy chain sequences.
#'
#' @param    db          data.frame with AIRR-format style columns.
#' @param    cell_id     column in \code{db} containing cell identifiers
#' @param    locus       column in \code{db} containing locus assignments
#' @param    sequence_id column in \code{db} containing locus assignments
#' @param    fields      Columns in \code{db}, in addition to \code{sample_id},
#'                       that should be used to group sequences to be
#'                       analyzed independently.
#'
#' @return   The input data.frame (\code{db}) with doublets removed.
#' @export
removeDoublets <- function(db, cell_id = "cell_id", locus = "locus",
                           sequence_id = "sequence_id", fields = NULL) {
    check <- alakazam::checkColumns(db,
                                    columns = c(cell_id, locus, sequence_id, fields))
    if (check != TRUE) { stop(check) }

    db_h <- db[isHeavyChain(db[[locus]]),, drop = FALSE] #IGH, TRB, TRD
    if (nrow(db_h) == 0) {
      message("`db` does not contain heavy chain sequences. Doublets check can't be performed.")
      return(db)
    }
    groups <- c(locus, fields)

    if (any(duplicated(db_h[[sequence_id]]))) {
        warning("Duplicated sequence ids in the input `db` (not considering `fields`).")
    }

    dup_ids <- sum(duplicated(db_h[, c(sequence_id, groups), drop = FALSE]))
    if (dup_ids > 0) {
        stop("Duplicated sequence ids found withing `c(locus,fields)` groups.")
    }

    doublets <- db_h %>%
                  group_by(!!!rlang::syms(c(cell_id, groups))) %>%
                  mutate(heavy_seqs_per_cell = n()) %>%
                  mutate(is_doublet = heavy_seqs_per_cell > 1) %>%
                  dplyr::filter(is_doublet) %>%
                  ungroup() %>%
                  select(!!!rlang::syms(c(fields, cell_id))) %>%
                  distinct()

    message("db size: ", nrow(db), " sequences.")
    message("number of heavy chain sequences: ", nrow(db_h))
    message("number of doublet cells: ", nrow(doublets))

    db %>% anti_join(doublets, by = c(fields, cell_id))
}

# Not used
# scQC <- function(db, cell_id="cell_id") {
#     db[["chain"]] <- getChain(db[["locus"]])
#     db[["tmp_scqc_row"]] <- 1:nrow(db)
#     db <- db %>%
#         group_by(sample_id, cell_id) %>%
#         mutate(paired_chain= all(c("VL", "VH") %in% chain)) %>%
#         ungroup() %>%
#         group_by(sample_id, chain, cell_id) %>%
#         mutate(num_sequences=n(),
#                num_isotypes=length(unique(c_call)),
#                cell_isotypes=paste(unique(c_call),sep=",",collapse=","),
#                max_count=max(consensus_count),
#                perc_of_cell_cons_count=100*consensus_count/sum(consensus_count),
#                is_most_abundant= consensus_count==max(consensus_count) & sum(consensus_count==max(consensus_count)) == 1,
#                cons_count_ratio=consensus_count/max_count) %>%
#         group_by(sample_id, chain, cell_id) %>%
#         do(cellQC(.)) %>%
#         ungroup()
#     qclog <- db %>%
#         arrange(tmp_scqc_row) %>%
#         select(
#                sample_id,
#                cell_id,
#                sequence_id,
#                paired_chain,
#                num_sequences,
#                num_isotypes,
#                cell_isotypes,
#                max_count,
#                perc_of_cell_cons_count,
#                is_most_abundant,
#                cons_count_ratio,
#                scqc_pass)
#     list(
#         "log"= qclog,
#         "pass"=sum(qclog$scqc_pass),
#         "fail"=sum(!qclog$scqc_pass)
#     )
# }
#
# cellQC <- function(df) {
#     if (length(unique(df$sample_id)) > 1) { stop ("Expectig data from one sample_id. Group by sample_id.")}
#     if (length(unique(df$chain)) > 1) { stop ("Expectig data from one chain. Group by chain.")}
#     df$scqc_pass <- FALSE
#     if (nrow(df)==1) {
#         df$scqc_pass <- T
#     } else {
#         df <- df %>%
#             arrange(desc(perc_of_cell_cons_count))
#         if (df$perc_of_cell_cons_count[1] >50) {
#             df$scqc_pass[1] <- TRUE
#         }
#         if (all(grepl("[LKlk]",df$chain)))  {
#             if (sum(df$perc_of_cell_cons_count[1:2]) > 75 & df$perc_of_cell_cons_count[2] > 25 ) {
#                 df$scqc_pass[1:2] <- TRUE
#             }
#         }
#     }
#     if (sum(df$scqc_pass) == 0) {
#         warning("No sequences selected for cell ", unique(df$cell_ide), "consensus_count: ", unique(df$consensus_count))
#     }
#     df %>%
#         arrange(tmp_scqc_row)
# }


#' findSingleCellDuplicates
#'
#' Finds sequences that share the same identifier and sequence between groups.
#'
#' @param    db    data.frame with AIRR-format style columns.
#' @param    fields      Columns in \code{db}, in addition to \code{sample_id},
#'                       that should be used to group sequences to be
#'                       analyzed independently.
#' @param    cell_id  column in \code{db} containing cell identifiers
#' @param    seq      column in \code{db} containing sequences to be compared
#' @param    sequence_id column in \code{db} containing sequence identifiers
#'
#' @return   A list with fields:
#'           \itemize{
#'             \item  \code{dups}:    a data.frame with the column \code{sc_duplicate}
#'                                    with values TRUE/FALSE to indicate whether the
#'                                    the row corresponds to a duplicated entry.
#'             \item  \code{fields}:  a data.frame showing the input fields used
#'             \item  \code{cell_id}: column in \code{db} containing cell identifiers
#'             \item  \code{seq}:     column in \code{db} containing sequence data
#'             \item  \code{sequence_id} column in \code{db} containin sequence identifiers
#'           }
#'
#' @export
findSingleCellDuplicates <- function(db, fields,
                                     cell_id = "cell_id",
                                     seq="sequence_alignment",
                                     sequence_id="sequence_id",
                                     mode=c("sequences", "cells")) {

    mode <- match.arg(mode)

    # Check for valid columns
    if (is.null(fields)) { stop("`fields` must be valid column name(s)")}
    columns <- c(fields,cell_id, seq, sequence_id)
    columns <- columns[!is.null(columns)]
    check <- alakazam::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }

    # Check that sequence_id are unique, because I will use this id
    # later to remove the duplicated sequences
    if (any(duplicated(db[[sequence_id]]))) {
        message("Duplicated `sequence_id` found. Expecting unique sequence identifiers. Using `row_id` internally.")
        db[["row_id"]] <- stri_join("row",1:nrow(db))
        sequence_id <- "row_id"
    }

    # Identify groups
    db[["group_idx"]] <- db %>%
                           dplyr::group_by(!!!rlang::syms(fields)) %>%
                           group_indices() %>%
                           paste0("gidx_", .)

    db <- db %>%
            select(!!!rlang::syms(c(fields, cell_id, seq, sequence_id,
                                    "group_idx")))

    # extract group names
    groups_table <- unique(db[, c("group_idx", fields), drop = FALSE])

    groups_table$group_name <- sapply(1:nrow(groups_table), function(g) {
        use <- colnames(groups_table) %in% c(fields) == T
        paste(as.matrix(groups_table[g, use, drop = FALSE]), collapse = "_")
    } )

    # Identify cells present in more than one group
    # that have same length sequences (to be able to compare them later
    # with .isDuplicateSeq (uses pairwiseDist))
    db_shared <- db %>%
                   mutate(seq_len = nchar(!!rlang::sym(seq))) %>%
                   select(!!!rlang::syms(c(columns, "group_idx", "seq_len"))) %>%
                   group_by(!!rlang::sym(cell_id)) %>%
                   mutate(is_shared_cell = length(unique(group_idx)) > 1) %>%
                   ungroup() %>%
                   dplyr::filter(is_shared_cell) %>%
                   select(-is_shared_cell)

    # Helper function to label as duplicate T/F those sequences
    # that have a 0 distance sequence in other group
    .isDuplicateSeq <- function(x, seq, sequence_id) {
        sequences <- x[[seq]]
        g <- x[["group_idx"]]
        # is_duplicate defaults to FALSE

        # Calculate distance between sequences
        # use gap=0 to treat gaps as Ns
        # and label as T those with a distance value of 0
        pwd <- pairwiseDist(sequences, dist_mat = getDNAMatrix(gap = 0)) == 0
        colnames(pwd) <- rownames(pwd) <- x[[sequence_id]]
        # don't use self comparison
        diag(pwd) <- NA
        pwd <- as.data.frame(pwd)
        pwd[[sequence_id]] <- rownames(pwd)

        dup_mat <- pivot_longer(pwd, -!!rlang::sym(sequence_id),
                                names_to = "group_seq") %>%
                    # filter(value == T) %>%
                    left_join(x[, c(sequence_id, "group_idx")],
                              by = sequence_id) %>%
                    rename(sequence_id_group = group_idx) %>%
                    left_join(x[, c(sequence_id, "group_idx")],
                              by = c("group_seq" = sequence_id)) %>%
                    dplyr::filter(sequence_id_group != group_idx) %>%
                    select(-group_seq, -sequence_id_group) %>%
                    group_by(!!!rlang::syms(c(sequence_id, "group_idx"))) %>%
                    summarize(n = sum(value, na.rm = T), .groups = "drop") %>%
                    ungroup() %>%
                    pivot_wider(id_cols = !!rlang::sym(sequence_id),
                                names_from = "group_idx", values_from = n)

        dup_mat[["is_duplicate_seq"]] <- as.logical(rowSums(dup_mat[, -1], na.rm = T) > 0)
        x %>% left_join(dup_mat, by = sequence_id )
    }

    # default to not duplicated
    db[["sc_duplicate"]] <- F
    duplicated_seq_id <- c()
    dups <- data.frame(matrix(ncol=length(fields), nrow=0, dimnames=list(NULL, fields)))

    if (nrow(db_shared) > 0 ) {
        db_shared <- db_shared %>%
            group_by(!!!rlang::syms(c(cell_id, "seq_len"))) %>%
            do(.isDuplicateSeq(., seq, sequence_id)) %>%
            ungroup() %>%
            select(-group_idx, -seq_len)
    }

    if ( sum(db_shared[["is_duplicate_seq"]],na.rm = T)>0 ) {

        if (mode == "sequences") {
            duplicated_seq_id <- db_shared %>%
                filter(is_duplicate_seq) %>%
                pull(!!rlang::sym(sequence_id))

            db[["sc_duplicate"]][db[[sequence_id]] %in% duplicated_seq_id] <- T

            dups <- db %>%
                left_join(db_shared %>% select(-is_duplicate_seq),
                          by = c(fields, sequence_id, seq, cell_id)) %>%
                rename(
                    sc_duplicate_group = group_idx
                )
        } else if (mode == "cells") {
            duplicated_cells <- db_shared %>%
                filter(is_duplicate_seq) %>%
                select(!any_of(c(seq, sequence_id, "is_duplicate_seq"))) %>%
                mutate(sc_duplicate = TRUE) %>%
                distinct()
            dups <- db %>%
                select(-sc_duplicate) %>%
                left_join(duplicated_cells,
                          by = c(fields, cell_id)) %>%
                rename(
                    sc_duplicate_group = group_idx
                )
            dups <- dups %>%
                mutate(sc_duplicate = replace_na(sc_duplicate, F))
        }

    }

    list(
        dups = dups,
        fields = groups_table,
        cell_id = cell_id,
        seq = seq,
        sequence_id = sequence_id
    )
}

# TODO: example
#' removeSingleCellDuplicates
#'
#' \code{removeSingleCellDuplicates} removes cells sharing the same identifier and sequence between groups.
#'
#' @param    db          data.frame with AIRR-format style columns.
#' @param    fields      Columns in \code{db}, in addition to \code{sample_id},
#'                       that should be used to group sequences to be
#'                       analyzed independently.
#' @param    cell_id  column in \code{db} containing cell identifiers
#' @param    seq      column in \code{db} containing sequences to be compared
#' @param    sequence_id column in \code{db} containing sequence identifiers
#' @param    mode     Use \code{sequences} to remove duplicated sequences and
#'                    \code{cells} to remove cells with duplicated sequences.
#' @return   A data.frame: a modified \code{db} without the duplicated sequences or cells.
#' @export
removeSingleCellDuplicates <- function(db, fields,
                                       cell_id = "cell_id",
                                       seq = "sequence_alignment",
                                       sequence_id = "sequence_id",
                                       mode = c("sequences", "cells" )) {

    mode <- match.arg(mode)

    check <- alakazam::checkColumns(db, columns = c(fields, cell_id, seq, sequence_id))
    if (check != TRUE) { stop(check) }

    # Check that sequence_id are unique, because I will use this id
    # later to remove the duplicated sequences
    if (any(duplicated(db[[sequence_id]]))) {
        message("Duplicated `sequence_id` found. Expecting unique sequence identifiers. Using `row_id` internally.")
        db[["row_id"]] <- stri_join("row",1:nrow(db))
        sequence_id <- "row_id"
    }

    dups <- findSingleCellDuplicates(db, fields, cell_id, seq, sequence_id)
    num_dups <- sum(dups[["dups"]][["sc_duplicate"]])

    if (num_dups > 0) {
        if (mode == "sequences") {
            message(num_dups, " sequences removed.")
            dups_true <- dups[["dups"]] %>% filter(sc_duplicate == T)
            db <- db %>% anti_join(dups_true[, sequence_id, drop = FALSE])
        } else {
            .is_sc_duplicate <- function(v){
                as.numeric(sum(v, na.rm = T) > 0)
            }
            # find cells
            dup_cells <- dups[["dups"]] %>%
                select(-!!rlang::sym(sequence_id), -!!rlang::sym(seq)) %>% # This is to count
                group_by(!!!rlang::syms(c(cell_id,fields))) %>%      # cells and not sequences
                summarize(across(starts_with("gidx_"), .is_sc_duplicate)) %>%
                ungroup() %>%
                mutate(num_samples = rowSums(select(., starts_with("gidx_")))) %>%
                filter(num_samples>0)
            message(nrow(dup_cells), " cells removed.")
            db <- db %>% anti_join(dup_cells[,c(fields, cell_id)])
        }
    } else {
        message("0 sequences/cells removed.")
    }

    if ("row_id" %in% colnames(db)) {
        db <- db %>% select(-row_id)
    }
    db
}

# TODO: document, return
# TODO: test
# dups <- findSingleCellDuplicates(db, groups="sample_id")
# counts <- singleCellSharingCounts(dups)
#' @export
singleCellSharingCounts <- function(dups) {
    if (nrow(dups[["dups"]]) < 1) {
      message("No duplicated cells found.")
      return(data.frame())
    }

    groups_table <- dups[["fields"]]
    groups <- setdiff(colnames(groups_table) , c("group_idx", "group_name"))
    cell <- dups[["cell_id"]]

    .getSmallerSize <- function(x, y, groups_sizes=groups_sizes_c) {
        sizes <- groups_sizes[c(x,y)]
        min(sizes, na.rm = T)
    }

    .rename <- function(gidxs) {
        name_map <- groups_table[["group_name"]]
        names(name_map) <- groups_table[["group_idx"]]
        name_map[gidxs]
    }

    # Find groups size. Duplicated cell ids in the same group
    # will be counted once
    groups_sizes <- dups[["dups"]] %>%
        select(!!!rlang::syms(c("sc_duplicate_group", cell))) %>%
        distinct() %>%
        group_by(sc_duplicate_group) %>%
        summarize(group_size = n()) %>% # group size is the number of unique cell ids
        ungroup()

    groups_sizes_c <- groups_sizes[["group_size"]]
    names(groups_sizes_c) <- groups_sizes[["sc_duplicate_group"]]

    # test <- data.frame(
    #  sample_id=c("s1", "s1", "s1","s2","s2", "s3"),
    #  cell_id  =c("c1","c1", "c2","c1","c3", "c3"),
    #  sequence_alignment = c("A", "T", "A", "A","A","T"),
    # sequence_id=c("seq1","seq2","seq3","seq4", "seq5", "seq6")
    # )
    # test_dups <- findSingleCellDuplicates(test, fields="sample_id",
    #                                      cell_id = "cell_id",
    #                                      seq="sequence_alignment",
    #                                      sequence_id="sequence_id")
    .is_sc_duplicate <- function(v) {
        as.numeric(sum(v, na.rm = T) > 0)
    }

    cell <- dups[["cell_id"]]
    sequence_id <- dups[["sequence_id"]]
    seq <- dups[["seq"]]
    counts <- dups[["dups"]] %>%
        select(-!!rlang::sym(sequence_id), -!!rlang::sym(seq)) %>% # This is to count
        group_by(sc_duplicate_group, !!rlang::sym(cell)) %>%      # cells and not sequences
        summarize(across(starts_with("gidx_"), .is_sc_duplicate)) %>%
        ungroup() %>%
        select(-!!rlang::sym(cell)) %>%
        group_by(sc_duplicate_group) %>%
        summarize(across(starts_with("gidx_"), sum)) %>%
        pivot_longer(-sc_duplicate_group, values_to="overlap_count")

    ## Add missing combinations, for a nice heatmap
    counts <- expand.grid(groups_table[["group_idx"]],
                          groups_table[["group_idx"]]) %>%
        rename(sc_duplicate_group = Var1,
               name = Var2) %>%
        left_join(counts) %>%
        mutate(overlap_count = ifelse(is.na(overlap_count), 0, overlap_count))

    ## Add self overlap, the diagonal
    counts <- counts %>%
        filter(sc_duplicate_group != name) %>%
        bind_rows(
            groups_sizes %>%
                mutate(name = sc_duplicate_group) %>%
                rename(overlap_count = group_size))

    counts <- counts %>%
        rowwise() %>%
        mutate(
            overlap_denominator = .getSmallerSize(sc_duplicate_group, name)
        ) %>%
        mutate(
            overlap_percent = 100*overlap_count/overlap_denominator,
            name = .rename(name),
            sc_duplicate_group = .rename(sc_duplicate_group)) %>%
        arrange(sc_duplicate_group, name)

    counts
}

# TODO: add caption
# TODO: document
# TODO: example
#' dups <- findSingleCellDuplicates(db, groups="sample_id")
#' counts <- singleCellSharingCounts(dups)
#' plotOverlapSingleCell(counts)
#' @export
plotOverlapSingleCell <- function(dup_count) {
    if (nrow(dup_count) < 1) {
      message("No duplicates found.")
      return(data.frame())
    }

    ggplot(dup_count, aes(sc_duplicate_group, name, group = name)) +
        geom_tile(aes(fill = overlap_percent)) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(x = "",y = "") +
        # annotate("segment", x = 0.5, xend = length(groups_table[["group_idx"]])+0.5,
        #          y = 0.5, yend = length(groups_table[["group_idx"]])+0.5, colour = "white") +
        geom_text(data = dup_count %>%
                           dplyr::filter(sc_duplicate_group == name),
                  aes(label = overlap_count),
                  angle = 0, color = "white", vjust = 0.5, hjust = 0.5,
                  na.rm = TRUE, size = 3)
}

##---------------------- TMP -------------------------- ##
# This function modified previous overlap heatmap function and is used to show, between two samples
# in single cell data, the overlap of cells which has same barcode and sequences
# find cells having same seq and cell barcode
#' singlecell_sharing_matrix(db)
#' @export
singlecell_sharing_matrix <- function(db,
                                      groups="sample_id",
                                      cell = "cell_id",
                                      seq="sequence_alignment",
                                      order_by=NULL) {
    # Check for valid columns
    if (is.null(groups)) { stop("`groups` must be a valid column name")}
    columns <- c(groups,cell, seq, order_by)
    columns <- columns[!is.null(columns)]
    check <- alakazam::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }

    # Identify groups
    db[["group_idx"]] <- db %>%
        dplyr::group_by(!!rlang::sym(groups)) %>%
        group_indices()

    # extract group names
    groups_table = unique(db[,c("group_idx",groups), drop = FALSE])

    groups_table$group_name <- sapply(1:nrow(groups_table), function(g) {
        use <- colnames(groups_table) %in% c(groups) == T
        paste(as.matrix(groups_table[g, use, drop = FALSE]), collapse = "_")
    } )

    num_groups = nrow(groups_table)

    # clonal assignment overlap matrices to be plotted
    overlap_matrix <- matrix(0, num_groups, num_groups)
    dimnames(overlap_matrix) <- list(groups_table$group_name, groups_table$group_name)

    overlap_count <- overlap_total <- overlap_percent <- overlap_matrix # matrix of overlap counts

    # calculate the number of overlapped cells
    for (ii in 1:nrow(groups_table)) {
        for (jj in 1:nrow(groups_table)) {

            ii_group_idx <- groups_table[["group_idx"]][ii]
            jj_group_idx <- groups_table[["group_idx"]][jj]

            # cell in groups i and j
            cells_ii = db[[cell]][db[["group_idx"]] == ii]
            cells_jj = db[[cell]][db[["group_idx"]] == jj]
            # seqs in groups i and j
            seqs_ii = db[[seq]][db[["group_idx"]] == ii]
            seqs_jj = db[[seq]][db[["group_idx"]] == jj]

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
    if (!is.null(order_by)) {
        overlap_count <- overlap_count[order_by,order_by]
        overlap_total <- overlap_total[order_by,order_by]
        overlap_percent <- overlap_percent[order_by,order_by]
    }
    return(createOverlap(overlap_count, overlap_total, overlap_percent))

}

# S4 class defining a Overlap object
setClass("Overlap",
         slots = c(overlap_count = "matrix",
                   overlap_total = "matrix",
                   overlap_percent = "matrix")
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
# plotOverlapSingleCell <- function(list_objOverlap="",
#                                   titleText="",
#                                   ROI=NULL,
#                                   order_by=NULL,
#                                   sizeDiag=3.5,
#                                   sizeCount=5,
#                                   sizePercent=4,
#                                   A=NULL,
#                                   B=NULL,
#                                   subject=""){
#
#     list_data_m <- list()
#
#     for(x in 1:length(list_objOverlap)){
#
#         objOverlap <- list_objOverlap[[x]]
#         if(!is.null(ROI)){
#             if(length(ROI)>1){
#                 roi <- match(ROI,row.names(objOverlap@overlap_count))
#                 roi <- roi[!is.na(roi)]
#             }else{
#                 roi <- which(grepl(ROI,row.names(objOverlap@overlap_count)))
#             }
#             objOverlap@overlap_count <- objOverlap@overlap_count[roi,roi]
#             objOverlap@overlap_total <- objOverlap@overlap_total[roi,roi]
#             objOverlap@overlap_percent <- objOverlap@overlap_percent[roi,roi]
#         }
#
#
#         data_m_count <- reshape::melt(objOverlap@overlap_count)
#         data_m_total <- reshape::melt(objOverlap@overlap_total)
#         data_m <- reshape::melt(objOverlap@overlap_percent)
#
#         #data_m$Y1 <- cut(data_m$value,breaks = c(0,0.00001,seq(0.00002,99.99998,by=1),99.99999,100),right = FALSE,include.lowest=TRUE)
#         data_m$Y1 <- cut(data_m$value,breaks = c(0,0.00001,1,2,3,4,6,10,20,30,99.99999,100),right = FALSE,include.lowest=TRUE)
#         #copa <- colorRampPalette(brewer.pal(5,"Greys"))(length(levels(data_m$Y1)))
#         copa <- colorRampPalette(brewer.pal(9,"Oranges"))(length(levels(data_m$Y1)))
#         names(copa) <- levels(data_m$Y1)
#
#         names(data_m)[3] <- "Percent"
#         data_m$Count <- data_m_count$value
#         data_m <- data_m[ is.finite(data_m$Percent), ]
#
#         #Format value to print
#         data_m$Value = round(data_m$Percent,2)
#         data_m$Value[ data_m$Value==0 ] = ""
#
#         # Non diagonal are counts/percent
#         #data_m$Value[ data_m$Value!="" ] = paste0( data_m$Count[ data_m$Value!="" ] , "\n(", data_m$Percent[ data_m$Value!="" ] , "%)")
#         data_m$Value[ data_m$Value!="" ] = data_m$Count[ data_m$Value!="" ]
#         data_m$Value[ as.vector(data_m$X1) == as.vector(data_m$X2) ] = ""
#
#         data_m$ValuePercent <- data_m$Value
#         data_m$ValuePercent[ data_m$Value!="" ] = data_m$Percent[ data_m$Value!="" ]
#
#         data_m$ValuePercent <-
#             unlist(
#                 sapply(data_m$ValuePercent, function(x){
#                     if(x!="") {
#                         return(paste0(format(round(as.numeric(x),1),nsmall=1),"%"))
#                     } else {
#                         return("")
#                     }
#                 }))
#
#         # Diagonals are counts
#         data_m$DiagonalCount = ""
#         data_m$DiagonalCount[as.vector(data_m$X1) == as.vector(data_m$X2)]= data_m$Count[as.vector(data_m$X1) == as.vector(data_m$X2)]
#
#         data_m$X1 <- factor(data_m$X1, levels=order_by)
#         data_m$X2 <- factor(data_m$X2, levels=order_by)
#
#         list_data_m[[x]] <- data_m
#     }
#
#
#     data_m1 <- list_data_m[[1]]
#     data_m2 <- list_data_m[[2]]
#
#     data_m <- data_m1
#     numbItems <- nrow(list_objOverlap[[1]]@overlap_count )
#     pos <- as.vector(upper.tri(matrix(NA, nrow=numbItems,ncol=numbItems,byrow=T)))
#
#
#     data_m$Value[pos] <-   data_m2$Value[pos]
#     data_m$ValuePercent[pos] <-   data_m2$ValuePercent[pos]
#     data_m$DiagonalCountLower <-   data_m$DiagonalCount
#     data_m$DiagonalCount <-   data_m2$DiagonalCount
#     data_m$Y1[pos] <-   data_m2$Y1[pos]
#     #sizeDiag=3
#     p <- ggplot(data=data_m, aes(X1, X2)) +
#         geom_tile(aes(fill = Y1), colour = "white") +
#         geom_text(aes(label = DiagonalCount), size=sizeDiag, angle = 50, color="white", vjust=-.5, hjust=0.5) +
#         geom_text(aes(label = DiagonalCountLower), size=sizeDiag, angle = 50, color="white", vjust=1.5, hjust=0.5) +
#         geom_text(aes(label = Value), size=sizeCount, color="black", vjust=0) +
#         #geom_text(aes(label = ValuePercent), size=sizePercent, color="black", vjust=1.4) +
#         scale_fill_manual(breaks=names(copa),values=copa, na.value="#FFFFFF", guide="none") +
#         theme(
#             axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black", size=10),
#             axis.text.y = element_text(colour = "black",size=10,hjust = 0),
#             axis.ticks.x=element_blank(),
#             axis.ticks.y=element_blank(),
#             text = element_text(size = 12, colour = "black"))+
#         annotate("segment", x = 0.5, xend = numbItems+.5, y = 0.5, yend = numbItems+.5, colour = "white") +
#         ggtitle(titleText) +
#         labs(x="",y="")+
#         theme(plot.title = element_text(lineheight=.8, face="bold"))+
#         #theme(plot.margin = unit(c(0,0,.25,.25), "cm")) +
#         theme(axis.title.y=element_text(vjust=1.3)) +
#         theme(axis.title.x=element_text(vjust=-.4))
#
#     return(p)
# }
