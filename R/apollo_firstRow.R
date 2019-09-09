#' Keeps only the first row for each individual
#' 
#' Given a multi-row input, keeps only the first row for each individual.
#' 
#' This a function to keep only the first row of an object per indidividual. It can handle multiple components, scalars, vectors and three-dimensional arrays (cubes).
#' The argument database MUST contain a column called 'apollo_sequence', which is created by \link{apollo_validateData}.
#' 
#' @param P List of vectors, matrices or 3-dim arrays. Likelihood of the model components (or other object).
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @return If \code{P} is a list, then it returns a list where each element has only the first row of each individual.
#'         If \code{P} is a single element, then it returns a single element with only the first row of each individual.
#'         The size of the element is changed only in the first dimension. If input is a scalar, then it returns a vector with the element repeated as many
#'         times as individuals in \code{database}. If the element is a vector, its length will be changed to the number of individuals. If the element is
#'         a matrix, then its first dimension will be changed to the number of individuals, while keeping the size of the second dimension. If the element 
#'         is a cube, then only the first dimension's length is changed, preserving the others.
#' @examples 
#' database <- data.frame(ID=rep(1:5, each=3), apollo_sequence=rep(1:3, 5))
#' apollo_inputs <- list(database=database)
#' attach(database)
#' 
#' P0 <- 0.5
#' apollo_firstRow(P0, apollo_inputs) 
#' 
#' P1 <- rep(c(0.1, 0.2, 0.3, 0.4, 0.5), each=3)
#' apollo_firstRow(P1, apollo_inputs) 
#' 
#' P2 <- matrix(rep(P1,10), nrow=15, ncol=10)
#' apollo_firstRow(P2, apollo_inputs) 
#' 
#' P3 <- array(rep(P1, 10*10), dim=c(15, 10, 10))
#' apollo_firstRow(P3, apollo_inputs)
#' 
#' P4 <- list(P0, P1, P2, P3)
#' apollo_firstRow(P4, apollo_inputs)
#' @export
apollo_firstRow=function(P, apollo_inputs){
  apollo_sequence <- apollo_inputs$database$apollo_sequence
  
  ### If P is a list
  if(is.list(P)){
    j=1
    while(j<=length(P)){
      isSca <- length(P[[j]])==1
      isVec <- is.vector(P[[j]]) && !isSca
      isMat <- is.matrix(P[[j]])
      isCub <- is.array(P[[j]]) && !isMat && length(dim(P[[j]]))==3
      if(isSca) P[[j]] = rep(P[[j]], sum(apollo_sequence==1))
      if(isVec|isMat) P[[j]] = subset(P[[j]],apollo_sequence==1)
      if(isCub) P[[j]] = P[[j]][(apollo_sequence==1),,,drop=FALSE]
      j=j+1
    }
    ### If P is a data.frame (besides being a list)
    if(is.data.frame(P)) P <- P[1:sum(apollo_sequence==1),]
  }else{ ### If P is not a list
    isSca <- length(P)==1
    isVec <- is.vector(P) && !isSca
    isMat <- is.matrix(P)
    isCub <- is.array(P) && !isMat && length(dim(P))==3
    if(isSca) P = rep(P, sum(apollo_sequence==1))
    if(isVec|isMat) P=subset(P,apollo_sequence==1)
    if(isCub) P=P[(apollo_sequence==1),,,drop=FALSE]
  }
  
  return(P)
}