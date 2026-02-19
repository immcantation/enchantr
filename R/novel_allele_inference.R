#' Create an Immcantation ...
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
novel_allele_inference_project <- function(path,...) {
    skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                              "templates", "project",
                              "novel_allele_inference_project_files")
    project_dir <- path
    if (!dir.exists(project_dir)) {
        message("Creating project_dir ", project_dir)
        dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    }
    project_files <- list.files(skeleton_dir, full.names = TRUE)
    file.copy(project_files, project_dir, recursive = TRUE)
}


#' Identify repeated SNP patterns in a sequence
#'
#' This function checks for specific repeated patterns of a newly identified SNP in a sequence.
#' If such a pattern is found, the sequence is masked and returned; otherwise, NA is returned.
#'
#' @param snp_string Character. SNP notation string (e.g., "123A>G").
#' @param germline_seq Character. The germline sequence to search for the pattern.
#' @return Character. Masked sequence if pattern is found, otherwise NA.
#' @examples
#' \dontrun{
#' Repeated_Read("123A>G", "ACTG...")
#' }
#' @export
Repeated_Read <- function(snp_string, germline_seq) {
   # Extract SNP position (NT), SNP nucleotide (SNP), and original nucleotide (OR_SNP)
   snp_pos <- as.numeric(gsub('([0-9]+).*', '\\1', snp_string))
   snp_nt <- gsub('.*>', '', snp_string)
   orig_nt <- gsub('[0-9]+([[:alpha:]]*).*', '\\1', snp_string)

   # Generate 4-mers around the SNP position (centered and shifted)
   kmer_seqs <- c(
      substr(germline_seq, snp_pos, snp_pos + 3),
      substr(germline_seq, snp_pos - 1, snp_pos + 2),
      substr(germline_seq, snp_pos - 2, snp_pos + 1),
      substr(germline_seq, snp_pos - 3, snp_pos)
   )

   # Build patterns for repeated SNPs and original nucleotide
   patterns <- c(
      paste0(rep(snp_nt, 3), orig_nt, collapse = ""),
      paste0(rep(snp_nt, 2), orig_nt, snp_nt, collapse = ""),
      paste0(snp_nt, orig_nt, rep(snp_nt, 2), collapse = ""),
      paste0(orig_nt, rep(snp_nt, 3), collapse = "")
   )
   pattern_regex <- paste0(patterns, collapse = '|')

   # Search for the pattern in the 4-mers
   match_idx <- grepl(pattern_regex, kmer_seqs)
   if (any(match_idx)) {
      # Mask the matched sequence: replace SNP with 'X', original with 'z'
      matched_seq <- kmer_seqs[match_idx]
      masked_seq <- gsub(snp_nt, 'X', gsub(orig_nt, 'z', matched_seq))
      return(masked_seq)
   } else {
      return(NA)
   }
}


#' Wrapper to check for repeated SNP patterns in a data.frame of novel alleles
#'
#' Applies the Repeated_Read function to each SNP substitution in each row of the novel alleles table.
#' Returns a logical vector indicating if any repeated SNP pattern was found for each row.
#'
#' @param novel_df Data.frame with columns 'nt_substitutions' (comma-separated SNPs) and 'germline_imgt' (sequence)
#' @return Logical vector: TRUE if a repeated SNP pattern is found in any substitution for the row, FALSE otherwise
#' \dontrun{
#' data("SampleNovel", package="tigger")
#' db <- SampleNovel
#' check_repeated_snp_patterns(db)
#' }
#' @export
check_repeated_snp_patterns <- function(novel_df) {
   vapply(seq_len(nrow(novel_df)), function(i) {
      snp_list <- strsplit(novel_df[['nt_substitutions']][i], ',')[[1]]
      germ_seq <- novel_df[['germline_imgt']][i]
      found <- vapply(snp_list, function(snp) {
         any(!is.na(Repeated_Read(snp, germ_seq)))
      }, logical(1))
      any(found)
   }, logical(1))
}


#' Calculate upper bounds for V gene germline ends (light chain)
#'
#' This function computes the 95th percentile (upper bound) of the germline end positions for each V gene,
#' and groups genes by their upper bound. Only V gene calls (not J) are considered. This is used for light chain
#' novel allele inference to set the region of interest for mutation analysis.
#'
#' @param data Data.frame with columns 'v_call' (V gene call) and 'v_germline_end' (end position of V germline alignment)
#' @return Named character vector: names are unique upper bounds, values are concatenated gene names for each bound
#' @details J gene sequences are not processed by this function. Only V gene calls without commas are included.
#' @examples
#' \dontrun{
#' data("ExampleDb", package="alakazam")
#' db <- ExampleDb
#' tigger_upper_bound(db)
#' }
#' @importFrom dplyr filter select mutate group_by summarise ungroup
#' @importFrom alakazam getGene
#' @export
tigger_upper_bound <- function(data) {
   # Filter to only V gene calls (no commas)
   v_data <- data %>%
      dplyr::filter(!grepl(',', .data[['v_call']])) %>%
      dplyr::select(v_call, v_germline_end)

   # Extract gene name
   v_data <- v_data %>%
      dplyr::mutate(GENE = alakazam::getGene(v_call, strip_d = FALSE))

   # Calculate 95th percentile (upper bound) for each gene
   v_data <- v_data %>%
      dplyr::group_by(GENE) %>%
      dplyr::mutate(RANGE = floor(quantile(v_germline_end, 0.95))) %>%
      dplyr::ungroup()

   # Group genes by their upper bound
   range_genes <- v_data %>%
      dplyr::select(GENE, RANGE) %>%
      dplyr::distinct() %>%
      dplyr::group_by(RANGE) %>%
      dplyr::summarise(GENES = paste0(GENE, "[*]", collapse = "|"), .groups = 'drop')

   # Return named vector: names are upper bounds, values are gene sets
   gene_range <- setNames(range_genes$GENES, range_genes$RANGE)
   return(gene_range)
}
