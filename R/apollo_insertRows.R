#' Inserts rows
#' 
#' Given a numeric object (scalar, vector, matrix or 3-dim array) inserts rows in the specified places.
#' 
#' In general, \code{r} should be longer than the number of rows in \code{v}, and \code{sum(r)=nrow(v)}.
#' If not, then a new object with as many rows as \code{r} will be returned. Old rows will be taken 
#' from \code{v} from the top down.
#' 
#' @param v Numeric scalar, vector, matrix or 3-dim array.
#' @param r Boolean vector. TRUE for inserting a row from \code{v}, \code{FALSE} to insert a new row with value \code{val}.
#' @param val Numeric scalar. Value that will fill new rows.
#' @return The same argument \code{v} but with the rows where \code{r==FALSE} removed.
#' @export
apollo_insertRows <- function(v, r, val){
  isVec <- is.vector(v)
  isMat <- is.matrix(v)
  isCub <- is.array(v) & length(dim(v))==3
  if(is.array(v)) nRows <- dim(v)[1] else nRows <- length(v)
  
  # If r is shorter than v
  if(length(r) < nRows){
    if(isVec) v <- v[1:sum(r)]
    if(isMat) v <- v[1:sum(r),,drop=FALSE]
    if(isCub) v <- v[1:sum(r),,,drop=FALSE]
  }
  
  if(isVec){
    v2    <- rep(val, length(r))
    v2[r] <- v
  }
  if(isMat){
    v2     <- matrix(val, nrow=length(r), ncol=ncol(v),
                     dimnames=list(NULL, colnames(v)))
    v2[r,] <- v
    
  }
  if(isCub){
    v2      <- array(val, dim=c(length(r), dim(v)[2], dim(v)[3]))
    dn <- dimnames(v)
    if(!is.null(dn)) dimnames(v2) <- list(NULL, dn[2], dn[3])
    v2[r,,] <- v
  }
  if(!isVec & !isMat & !isCub) v2 <- v
  
  return(v2)
}