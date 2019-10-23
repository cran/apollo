#' Writes an entry to \code{apolloLog}
#'
#' Writes an entry to the apolloLog, which lives inside \code{apollo_inputs}.
#'
#' The variable \code{apolloLog} is an environment created inside \code{apollo_inputs}
#' by \code{apollo_validateInputs}, but re-set by apollo_estimate.
#' As an environment, it can be modified in place, i.e. all changes done within this 
#' function are recorded in apolloLog, even if it belongs to another environment.
#' @param title Character. Title of the log entry.
#' @param content Content of the log entry. Can be a single element or a list.
#'                Each element will be converted to character using \code{print},
#'                and concatenated with a line feed in between.
#' @param apolloLog Environment. It contains the character vectors of titles and content.
#' @importFrom utils capture.output
#' @return TRUE if writing was succesful, FALSE if not.
apollo_addLog <- function(title="", content="", apolloLog){
  ### Get apolloLog from apollo_inputs
  #apolloLog <- apollo_getLog()
  #if(length(apolloLog)==1 && is.na(apolloLog)) return(FALSE)
  
  ### Validation
  if(!is.environment(apolloLog)) return(FALSE)
  
  ### Convert title and content to text if necessary
  # apollo_addLog(title="lista", content=list(matrix(1:9,3,3), "hi", "mate", ""))
  # x <- apollo_getLog()
  # cat(apollo_getLog(TRUE))
  if(length(title)!=1 || !is.character(title)){
    title <- capture.output(print(title))
  }
  if(is.list(content)){
    content[which(content %in% "")] <- NULL
    tmp <- ""
    for(i in content){
      if(is.character(i)){
        tmp <- paste0(tmp, ifelse(nchar(tmp)>0,"\n",""), paste(i, collapse="\n"))
      } else{
        x <- paste(capture.output(print(i)), collapse="\n")
        tmp <- paste0(tmp, ifelse(nchar(tmp)>0,"\n",""), x, "\n")
      }
      #tmp <- paste0(tmp, "\n")
    }
    # Remove las character if its "\n"
    lastChar <- substring(tmp, nchar(tmp), nchar(tmp))
    if(lastChar=="\n") content <- substring(tmp, 1, nchar(tmp)-1) else content <- tmp
  } else{
    if(is.character(content)){
      content <- paste0(paste(content, collapse="\n"))
    }
    if(!is.character(content)){
      content <- paste0(paste(capture.output(print(content)), collapse="\n"))
    }
  }
  
  ### Check that new input is different to previous ones
  if(title %in% apolloLog$title){
    if(content %in% apolloLog$content) return(FALSE)
  }
  
  
  ### Add new entry to apolloLog
  apolloLog$title <- c(apolloLog$title, title)
  apolloLog$content <- c(apolloLog$content, content)
  
  return(TRUE)
}