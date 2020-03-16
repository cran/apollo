#' Returns the log of Apollo
#'
#' Returns the apolloLog variable either as a list or as text.
#'
#' The variable apolloLog is a list whose elements are character vectors with two elements.
#' The first element is the title of the entry, and the second element is content of the entry.
#' ApolloLog lives in the namespace environment of the Apollo package.
#' @param apolloLog Environment. It contains the character vectors of titles and content.
#' @param book Positive scalar integer. Book number inside the log (default is 1).
#' @return A list or a scalar character variable.
apollo_printLog <- function(apolloLog, book=1){
  ### Check that apolloLog is not corrupted
  test <- !is.environment(apolloLog)
  test <- test || !exists("title", envir=apolloLog)
  test <- test || !exists("content", envir=apolloLog)
  test <- test || !is.list(apolloLog$title)
  test <- test || !is.list(apolloLog$content)
  test <- test || length(apolloLog$title)!=length(apolloLog$content)
  test <- test || length(apolloLog$title) < book
  if(test) return("")
  
  ### Convert to text
  tmp <- c()
  n <- min( length(apolloLog$title[[book]]), 
            length(apolloLog$content[[book]]) )
  if(n>0) for(i in 1:n){
    tmp <- paste0(tmp, 
                  apolloLog$title[[book]][i], "\n", 
                  apolloLog$content[[book]][i], "\n",
                  ifelse(i<n, "\n", ""))
  } 
  return(tmp)
}