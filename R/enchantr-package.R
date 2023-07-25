#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' enchantr
#'
#' Reports
#'
#' For additional details regarding the use of the \code{enchantr} package see the
#' vignettes:\cr
#' \code{browseVignettes("enchantr")}
#'
#'
#' @section  Reports:
#' \itemize{
#'   \item  \link{chimera_analysis_project}:    Create an Immcantation Chimera detection project
#'   \item  \link{collapse_duplicates_project}: Create an Immcantation Collapse duplicates project
#'   \item  \link{contamination_project}:       Create an Immcantation Detect contamination project
#'   \item  \link{find_threshold_project}:      Create an Find threshold project
#'   \item  \link{define_clones_project}:       Create an Immcantation Define clones project
#'   \item  \link{dowser_lineage_project}:      Create an Immcantation Dowser project
#'   \item  \link{single_cell_qc_project}:      Create an Immcantation sc QC project
#' }
#'
#'
#' @name     enchantr
#' @docType  package

#' @import      alakazam
#' @importFrom  ComplexHeatmap make_comb_mat UpSet
#' @importFrom  foreach foreach %dopar% registerDoSEQ
#' @importFrom  doParallel  registerDoParallel
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
#'                      make_line_graph vcount as_long_data_frame
#' @importFrom  iterators   icount
#' @importFrom  plotly  ggplotly
#' @importFrom  RColorBrewer brewer.pal
#' @importFrom  rlang       sym syms
#' @importFrom  stringi stri_trim_both stri_split_fixed stri_join
#' @importFrom  stringr str_replace str_replace_all
#' @importFrom  tidyr pivot_longer pivot_wider
#' @importFrom  xfun       in_dir
#' 
## usethis namespace: end
NULL
