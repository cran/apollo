#' Keeps only some rows
#' 
#' Given a numeric object (scalar, vector, matrix or 3-dim array) keeps only the specified rows.
#' 
#' @param v Numeric scalar, vector, matrix or 3-dim array.
#' @param r Boolean vector. As many elements as rows in \code{v}. TRUE for keeping the row. \code{FALSE} to drop it.
#' @return The same argument \code{v} but with the rows where \code{r==FALSE} removed.
#' @export
apollo_keepRows <- function(v, r){
  # If v is scalar, but r isn't, don't change anything
  if(length(v)==1){
    if(length(r)==1) return(v[r]) else return(v)
  }
  
  if(is.data.frame(v)) return(v[r,])
  if(is.array(v) && length(dim(v))==3) return(v[r,,,drop=FALSE])
  if(is.matrix(v)) return(v[r,,drop=FALSE])
  if(!is.array(v) && is.vector(v)) return(v[r])
  
  return(v)
}