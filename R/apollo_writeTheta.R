#' Writes the vector [beta,ll] to a file called \code{modelname_iterations.csv}
#' @param beta vector of parameters to be written (including fixed ones).
#' @param ll scalar representing the log-likelihood of the whole model.
#' @param modelName Character. Name of the model.
#' @param scaling Numeric vector of scales applied to beta
#' @param outDir Scalar character. Name of output directory
#' @param apollo_beta Named numeric vector of starting values.
#' @return Nothing.
#' @export
apollo_writeTheta <- function(beta, ll, modelName, scaling=NULL, 
                              outDir=NULL, apollo_beta=NULL){
  
  ### Scale parameters
  test <- is.numeric(scaling) && is.vector(scaling) && !is.null(names(scaling))
  if(test) beta[names(scaling)] <- beta[names(scaling)]*scaling
  
  ### Set output directory
  test <- !is.null(outDir) && is.character(outDir) && length(outDir)==1
  if(!test) outDir <- paste0(getwd(),'/')
  
  ### Sort parameters in the same order than apollo_beta
  bNam <- names(apollo_beta)
  test <- is.character(bNam) && all(bNam %in% names(beta))
  if(test) beta <- beta[bNam]
  
  ### Initialise
  fileName <- paste0(outDir,modelName,"_iterations.csv")
  tmp <- matrix(c(beta,ll), nrow=1)
  
  if(file.exists(fileName)){
    # If file already exists, append
    tryCatch(utils::write.table(tmp, file=fileName, append=TRUE, sep=",", 
                                col.names=FALSE, row.names=FALSE),
             error=function(e) cat("Current iteration could not be written to ",
                                   fileName, ".\n", sep="") )
  } else {
    # If file does not exist, write afresh with headings
    colnames(tmp) <- c(names(beta),'logLike')
    tryCatch(utils::write.table(tmp, file=fileName, sep=',', row.names=FALSE, 
                                append=FALSE),
             error=function(e) cat("Initial iteration could not be written to ",
                                   fileName, ".\n", sep="") )
  }
}
