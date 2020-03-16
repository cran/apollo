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
#' @param book Positive scalar integer. Book number inside the log (default is 1).
#' @importFrom utils capture.output
#' @return TRUE if writing was succesful, FALSE if not.
apollo_addLog <- function(title="", content="", apolloLog, book=1){
  ### Get apolloLog from apollo_inputs
  #apolloLog <- apollo_getLog()
  #if(length(apolloLog)==1 && is.na(apolloLog)) return(FALSE)
  
  ### Validation
  if(!is.environment(apolloLog)) return(FALSE)
  
  ### Convert title and content to text if necessary
  if(length(title)!=1 || !is.character(title)){
    title <- capture.output(print(title))
  }
  if(is.list(content)){ # if content is a list
    content[which(content %in% "")] <- NULL
    tmp <- ""
    for(i in content){
      if(is.character(i)){
        tmp <- paste0(tmp, ifelse(nchar(tmp)>0,"\n",""), paste(i, collapse="\n"))
      } else{
        x <- paste(capture.output(print(i)), collapse="\n")
        tmp <- paste0(tmp, ifelse(nchar(tmp)>0,"\n",""), x, "\n")
      }
    }
    # Remove last character if its "\n"
    lastChar <- substring(tmp, nchar(tmp), nchar(tmp))
    if(lastChar=="\n") content <- substring(tmp, 1, nchar(tmp)-1) else content <- tmp
  } else{ # if content is not a list
    if(is.character(content)){
      content <- paste0(paste(content, collapse="\n"))
    }
    if(!is.character(content)){
      content <- paste0(paste(capture.output(print(content)), collapse="\n"))
    }
  }
  
  ### Create title and content list if necessary
  test <- !exists("title", envir=apolloLog) || !is.list(apolloLog$title)
  if(test) apolloLog$title <- list()
  #if(length(apolloLog$title)<book) apolloLog$title[[book]] <- ""
  test <- !exists("content", envir=apolloLog) || !is.list(apolloLog$content)
  if(test) apolloLog$content <- list()
  #if(length(apolloLog$content)<book) apolloLog$content[[book]] <- ""
  
  ### Check that new input is different to previous ones
  test <- length(apolloLog$title)>=book
  test <- test && length(apolloLog$content)>=book
  test <- test && ( title   %in% apolloLog$title[[book]] )
  test <- test && ( content %in% apolloLog$content[[book]] )
  if(test) return(FALSE)
  
  ### Add new entry to apolloLog
  test <- length(apolloLog$title)>=book
  test <- test && length(apolloLog$content)>=book
  if(!test){ # create book
    apolloLog$title[[book]] <- title
    apolloLog$content[[book]] <- content
  } else { # add to existing book
    apolloLog$title[[book]] <- c(apolloLog$title[[book]], title)
    apolloLog$content[[book]] <- c(apolloLog$content[[book]], content)
  }
  
  
  return(TRUE)
}