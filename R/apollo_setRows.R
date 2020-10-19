#' Sets specified rows to a given value
#' 
#' Given a numeric object (scalar, vector, matrix or 3-dim array) sets a subset of rows to a given value.
#' 
#' @param v Numeric scalar, vector, matrix or 3-dim array. Rows of this object will be replaced by \code{val} and 
#' @param r Boolean vector. As many elements as rows in \code{v}. TRUE for replacing that row, FALSE for not changing it.
#' @param val Numeric scalar. Value to which the specified rows must be set to.
#' @return The same argument \code{v} but with the rows where \code{r==TRUE} set to \code{val}.
#' @export
apollo_setRows <- function(v, r, val){
  # If v is scalar, but r isn't, don't change anything
  if(length(v)==1){
    if(length(r)==1 && r) return(val) else return(v)
  }
  
  if(is.array(v) && length(dim(v))==3) v[r,,] <- val
  if(is.matrix(v)) v[r,] <- val
  if(!is.array(v) && is.vector(v)) v[r] <- val
  
  return(v)
}