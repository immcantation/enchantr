#' Create an Immcantation ...
#'
#' From RStudio, use the New Project wizard: File > New Project >
#' New Directory > then select  Immcantation ...
#' to create the skeleton of an Immcantation ... project
#' @param  path path to the directory where the project will be created
repertoire_analysis_project <- function(path,...) {
  skeleton_dir <- file.path(system.file(package = "enchantr"), "rstudio",
                            "templates", "project",
                            "repertoire_analysis_project_files")
  project_dir <- path
  if (!dir.exists(project_dir)) {
    message("Creating project_dir ", project_dir)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_files <- list.files(skeleton_dir, full.names = TRUE)
  file.copy(project_files, project_dir, recursive = TRUE)
}


#' Plot heatmap of clonal overlap
#'
#' Plot a matrix to visualize the the number of clones shared between groups.
#'
#' @param db              Changeo db data.frame
#' @param group           Vector with column names for grouping. Overlap will
#'                        be calculated across the groups. e.g \code{SAMPLE}
#' @param heatmap_colors  Vector of colors representing low and high values on
#'                        the heatmap and the diagonal. Default is
#'                        c("white","orange", "grey80")
#' @param xlab            Text to be used as x axis title
#' @param ylab            Text to be used as y axis title
#' @param geom_text_size  Plot text size
#' @export
#' 
ClonalOverlap_plot <- function(db, group="sample_id",
                               heatmap_colors=c("white","orange", "grey80"),
                               xlab=NULL, ylab=NULL,
                               geom_text_size=3){
  
  if (length(unique(db[[group]])) == 1) {
    warning("One group only. Can't look for overlaps.")
    return(NULL)
  }
  
  unique_group <-str_sort(unique(db[[group]]),numeric=TRUE)
  
  sample_x <- c()
  sample_y <- c()
  nclones_x <- c()
  nclones_y <- c()
  n_shared_clones <- c()
  perc_shared_y_vec <- c()
  
  for (x in unique_group){
    for (y in unique_group){
      sample_x <- c(sample_x,x)
      sample_y <- c(sample_y,y)
      clones_x = unique(db[["clone_id"]][db[[group]] == x])
      clones_y = unique(db[["clone_id"]][db[[group]] == y])
      n_clone_x = length(clones_x)
      n_clone_y = length(clones_y)
      num_shared_clones <- length(intersect(clones_x, clones_y))
      perc_shared_y <- round ( num_shared_clones / n_clone_y * 100 , 1)
      nclones_x <- c(nclones_x, n_clone_x)
      nclones_y <- c(nclones_y, n_clone_y)
      n_shared_clones <- c(n_shared_clones, num_shared_clones)
      perc_shared_y_vec <- c(perc_shared_y_vec,perc_shared_y )
    }
  }
  
  df <- data.frame(sample_x, sample_y, nclones_x, nclones_y, n_shared_clones, perc_shared_y_vec)
  colnames(df) <- c("sample_x", "sample_y","Sample_x_Clone_Numbers", "Sample_y_Clone_Numbers", "Overlap_Numbers", "Overlap_Percentage")
  
  df$Overlap_Percentage <- as.numeric(df$Overlap_Percentage)
  
  df <- df %>% mutate(overlap=paste0(Overlap_Numbers, "\n(", Overlap_Percentage, "%)"))
  
  df <- df %>%
    mutate(
      Overlap_Percentage = ifelse(sample_x == sample_y,
                                  NA_real_,  # keeps column numeric
                                  Overlap_Percentage)
    )
  
  df <- df %>%
    mutate(overlap = ifelse(Overlap_Numbers == 0, "", overlap))
  
  p <- ggplot(df, aes(x=sample_x, y=sample_y)) +
    geom_tile(aes(fill = Overlap_Percentage, text = overlap)) +
    theme_enchantr()  + 
    scale_fill_gradient(name = "Overlap Percentage", 
                        low = heatmap_colors[1], high = heatmap_colors[2], na.value = heatmap_colors[3]) + 
    xlab(xlab) + ylab(ylab) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (geom_text_size > 0) {
    p <- p + geom_text(aes(label = overlap), size = geom_text_size, na.rm = TRUE)
  }
  
  return(p)
}


#' Plot heatmap of bootstrapping clonal overlap 
#'
#' Plot a matrix to visualize the the average clonal overlap percentage in bootstrapping draws between groups.
#'
#' @param bootstrap_db    bootstrapping data.frame after calling estimateAbundance 
#' @param group           Vector with column names for grouping. Overlap will
#'                        be calculated across the groups. e.g \code{SAMPLE}
#' @param heatmap_colors  Vector of colors representing low and high values on
#'                        the heatmap and the diagonal. Default is
#'                        c("white","orange", "grey80")
#' @param xlab            Text to be used as x axis title
#' @param ylab            Text to be used as y axis title
#' @param geom_text_size  Plot text size
#' @export
#' 
Bootstrapping_ClonalOverlap_plot <- function(bootstrap_db, group="sample_id",
                               heatmap_colors=c("white","orange", "grey80"),
                               xlab=NULL, ylab=NULL,
                               geom_text_size=3){
  
  if (length(unique(bootstrap_db[[group]])) == 1) {
    warning("One group only. Can't look for overlaps.")
    return(NULL)
  }
  
  unique_group <-str_sort(unique(bootstrap_db[[group]]),numeric=TRUE)
  
  sample_x_vector <- c()
  sample_y_vector <- c()
  average_overlap_y_vector <- c()
  
  n<-length(unique_group)
  
  # Pairwise clonal overlap calculation
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      x <- unique_group[i]
      y <- unique_group[j]
      sample_x_vector <- c(sample_x_vector,x)
      sample_y_vector <- c(sample_y_vector,y)
      res <- bootstrap_shared_percentage(bootstrap_db%>%filter(.data[[group]]%in%c(x,y)), x, y, group)
      average_overlap_x <-round(res$average_overlap_x*100,2)
      average_overlap_y <- round(res$average_overlap_y*100,2)
      average_overlap_y_vector <- c(average_overlap_y_vector, average_overlap_y)
      # fill the cell with flipped x and y value
      sample_x_vector <- c(sample_x_vector,y)
      sample_y_vector <- c(sample_y_vector,x)
      average_overlap_y_vector <- c(average_overlap_y_vector, average_overlap_x)
    }
  }
  
  df <- data.frame(sample_x_vector, sample_y_vector, average_overlap_y_vector)
  colnames(df) <- c("sample_x", "sample_y", "Average_Overlap_Percentage")
  df <- df %>% mutate(overlap=paste0(Average_Overlap_Percentage, "%"))
  
  p <- ggplot(df, aes(x=sample_x, y=sample_y)) +
    geom_tile(aes(fill = Average_Overlap_Percentage, text = overlap)) +
    theme_enchantr()  + 
    scale_fill_gradient(name = "Overlap Percentage", 
                        low = heatmap_colors[1], high = heatmap_colors[2], na.value =           heatmap_colors[3]) + 
    xlab(xlab) + ylab(ylab) + 
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (geom_text_size > 0) {
    p <- p + geom_text(aes(label = overlap), size = geom_text_size, na.rm = TRUE)
  }
  return(p)
}


#' Calculate clonal overlap between two samples from bootstrapping results
#'
#' @param bootstrap_db    bootstrapping data.frame after calling estimateAbundance 
#' @param x               sample x
#' @param y               sample y
#' @param group           Vector with column names for grouping. Overlap will
#'                        be calculated across the groups. e.g \code{SAMPLE}
#' 
bootstrap_shared_percentage<- function(bootstrap_db, x, y, group='sample_id'){
  # Select bootstrapping table of only sample x and sample y and seen clones. 
  compare_db <- bootstrap_db %>% 
    filter(.data[[group]]%in%c(x,y)) %>% 
    filter(! str_starts(clone_id, 'U'))  # filter out inferred unseen clones
  
  seen_clones_x <- compare_db %>% filter(.data[[group]]==x) %>% pull(clone_id)
  seen_clones_y <- compare_db %>% filter(.data[[group]]==y) %>% pull(clone_id)
  shared_clones <- intersect(seen_clones_x, seen_clones_y)
  
  if(length(shared_clones)==0){
    result <- list(
      average_overlap_x = 0,
      average_overlap_y = 0
    )
    return(result)
  }
  else{
    compare_db_x <- compare_db %>% 
      filter(clone_id %in% shared_clones) %>% 
      arrange(clone_id) %>%           # order by clone_id
      filter(.data[[group]]==x) %>%
      select(3:ncol(compare_db))%>%
      replace(.>0,1)                  # replace the value greater than 0 to 1
    
    compare_db_y <- compare_db %>% 
      filter(clone_id %in% shared_clones) %>% 
      arrange(clone_id) %>%
      filter(.data[[group]]==y) %>%
      select(3:ncol(compare_db))%>%
      replace(.>0,1)
    
    new_db <- compare_db_x & compare_db_y  # logical & operation on the two bootstrapping table from sample x and sample y 
    
    new_db <- new_db*1       # convert logic to integer
    
    overlap_vec <- colSums(new_db)   # vector of overlap count per bootstrapping draw
    
    # vector of clone observation count per bootstrapping draw
    clone_observation_x <- colSums(compare_db %>% 
                                     filter(.data[[group]]==x) %>%
                                     select(3:ncol(compare_db))%>%
                                     replace(.>0,1))
    
    clone_observation_y <- colSums(compare_db %>% 
                                     filter(.data[[group]]==y) %>%
                                     select(3:ncol(compare_db))%>%
                                     replace(.>0,1))
    
    # Sum overlap percentage per bootstrapping draw and average it by number of bootstrapping draws
    average_overlap_x<-sum(dplyr::if_else(clone_observation_x == 0,
                                          0,
                                          overlap_vec / clone_observation_x)) /(ncol(bootstrap_db) - 2)
    average_overlap_y<-sum(dplyr::if_else(clone_observation_y == 0,
                                          0,
                                          overlap_vec / clone_observation_y)) /(ncol(bootstrap_db) - 2)
    
    # average number of overlap 
    result <- list(
      average_overlap_x       = average_overlap_x,
      average_overlap_y = average_overlap_y
    )
    
    return(result)
  }
  
}