#' Keeps only the first row for each individual
#' 
#' Given a multi-row input, keeps only the first row for each individual.
#' 
#' This a function to keep only the first row of an object per indidividual. It can handle multiple types of components, including scalars, vectors and three-dimensional arrays (cubes).
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
#' @export
apollo_firstRow=function(P, apollo_inputs){
  ### Fetch functionality
  functionality <- tryCatch(get('functionality', envir=parent.frame(), inherits=FALSE),
                            error=function(e) NULL)
  if(!is.null(functionality) && functionality %in% c('preprocess', 'report')) return(P)
  
  apollo_sequence <- apollo_inputs$database$apollo_sequence
  
  if(is.data.frame(P)){
    if(nrow(P)==length(apollo_sequence)){
      return(subset(P,apollo_sequence==1))
    }else{
      stop("INPUT ISSUE - The data.frame passed to apollo_firstRow does not have the same number of rows as the database!")
    } 
  }
  
  ### If P is a list
  if(is.list(P)){
    for(j in 1:length(P)){
      isSca <- length(P[[j]])==1
      isVec <- is.vector(P[[j]]) && !isSca
      isMat <- is.matrix(P[[j]])
      isCub <- is.array(P[[j]]) && !isMat && length(dim(P[[j]]))==3
      condition <- isSca || (isVec && length(P[[j]])==length(apollo_sequence)) || (isMat && nrow(P[[j]])==length(apollo_sequence)) || (isCub && nrow(P[[j]])==length(apollo_sequence))
      if(!condition) stop("INPUT ISSUE - The object passed to apollo_firstRow does not have the same number of entries as the database!")
      if(isSca) P[[j]] = rep(P[[j]], sum(apollo_sequence==1))
      if(isVec|isMat) P[[j]] = subset(P[[j]],apollo_sequence==1)
      if(isCub) P[[j]] = P[[j]][(apollo_sequence==1),,,drop=FALSE]
    }
  }else{ ### If P is not a list
    isSca <- length(P)==1
    isVec <- is.vector(P) && !isSca
    isMat <- is.matrix(P)
    isCub <- is.array(P) && !isMat && length(dim(P))==3
    condition <- isSca || (isVec && length(P)==length(apollo_sequence)) || (isMat && nrow(P)==length(apollo_sequence)) || (isCub && nrow(P)==length(apollo_sequence))
    if(!condition) stop("INPUT ISSUE - The object passed to apollo_firstRow does not have the same number of entries as the database!")
    if(isSca) P = rep(P, sum(apollo_sequence==1))
    if(isVec|isMat) P=subset(P,apollo_sequence==1)
    if(isCub) P=P[(apollo_sequence==1),,,drop=FALSE]
  }
  
  return(P)
}