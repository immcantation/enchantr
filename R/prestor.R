# Borrowed from prestoR
# https://bitbucket.org/kleinstein/prestor/raw/8844ac4837a51fadbb91bdcb3a1debf033ec7403/R/Core.R

#' Parse console output from a pRESTO pipeline
#'
#' @param    log_file    console output filename
#' 
#' @return   data.frame  with columns (step, task, field, value)
#' 
#' @export
loadConsoleLog <- function(log_file) {
    ## DEBUG
    # log_file <- file.path(IN_PATH, n, "Pipeline.log")
    
    # Parse console output
    log_text <- scan(log_file, what=character(), sep="\n", quiet=TRUE)
    log_text <- stri_trim_both(log_text)
    log_text <- log_text[!grepl("(^PROGRESS>)|(^END>)", log_text)]
    log_df <- as.data.frame(stri_split_fixed(log_text, "> ", 2, simplify=TRUE), stringsAsFactors=FALSE)
    names(log_df) <- c("field", "value")
    
    # Assign steps and tasks
    log_df$step <- NA
    log_df$task <- NA
    x <- c(which(log_df[1] == "START"), nrow(log_df))
    for (i in 2:length(x)) {
        r <- x[i-1]:x[i]
        log_df$step[r] <- i - 1
        if ("COMMAND" %in% log_df$field[r]) {
            n <- paste(log_df$value[x[i-1]], 
                       log_df$value[r][log_df$field[r] == "COMMAND"],
                       sep="-")
        } else {
            n <- log_df$value[x[i-1]]
        }
        log_df$task[r] <- n
    }
    
    # Remove START field
    #log_list <- dlply(log_df, .(step))
    log_df <- subset(log_df, !(field %in% c("START", "COMMAND")), select=c(step, task, field, value))
    log_df <- transform(log_df, step=as.numeric(step), task=factor(task, levels=unique(task)))
    
    return(log_df)
}

