#' Returns the log of Apollo
#'
#' Returns the apolloLog variable either as a list or as text.
#'
#' The variable apolloLog is a list whose elements are character vectors with two elements.
#' The first element is the title of the entry, and the second element is content of the entry.
#' ApolloLog lives in the namespace environment of the Apollo package.
#' @param apolloLog Environment. It contains the character vectors of titles and content.
#' @return A list or a scalar character variable.
apollo_printLog <- function(apolloLog){
  ### Get apolloLog from parent environments, hopefully
    # from the execution environment of apollo_estimate
  #apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE )$apolloLog, 
  #                      error=function(e){
  #                        warning("Could not retrieve apolloLog")
  #                        return(NA)})
  
  if((length(apolloLog)==1 && is.na(apolloLog)) || is.null(apolloLog)) return(NA)
  if(!is.environment(apolloLog) || 
     (is.environment(apolloLog) && length(apolloLog$title)!=length(apolloLog$content))){
    warning("apolloLog is corrupted")
    return(NA)
  }
  if(!exists("title", envir=apolloLog) || !exists("content", envir=apolloLog)){
    apolloLog$title <- c()
    apolloLog$content <- c()
  }
  
  ### Convert to text
  n <- min(length(apolloLog$title), length(apolloLog$content))
  if(n>0){
    tmp <- ""
    for(i in 1:n) tmp <- paste0(tmp, 
                                apolloLog$title[i], "\n", 
                                apolloLog$content[i], "\n",
                                ifelse(i<n, "\n", ""))
    return(tmp)
  } else return("Apollo log is empty\n")
}