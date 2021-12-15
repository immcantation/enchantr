#' enchantr
#' 
#' Reports
#' 
#' @name     enchantr
#' @docType  package

#' @import      alakazam
#' @importFrom  foreach foreach %dopar% registerDoSEQ
#' @import      dplyr
#' @import      data.table
#' @import      DT
#' @import      ggplot2
#' @import      ggraph
#' @importFrom  airr    read_rearrangement write_rearrangement
#' @importFrom  igraph  E V V<- E<- 
#'                      add_edges add_vertices
#'                      as_edgelist  
#'                      decompose degree
#'                      graph_from_data_frame graph_from_edgelist
#'                      make_line_graph vcount
#' @importFrom  RColorBrewer brewer.pal
#' @importFrom  stringi stri_trim_both stri_split_fixed
#' @importFrom  stringr str_replace str_replace_all
#' @importFrom  tidyr pivot_longer pivot_wider
NULL