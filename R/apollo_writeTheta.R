#' Writes the vector [beta,ll] to a file called modelname_iterations.csv
#' @param beta vector of parameters to be written.
#' @param ll scalar representing the loglikelihood of the whole model.
#' @param modelName Character. Name of the model.
#' @return Nothing.
#' @export
apollo_writeTheta <- function(beta, ll, modelName){
  fileName <- paste(modelName,'iterations.csv',sep='_')
  tmp <- matrix(c(beta,ll ),nrow=1)
  
  if(file.exists(fileName)){
    tryCatch( utils::write.table(tmp,file=fileName, append=TRUE, sep=',', col.names=FALSE, row.names=FALSE),
              error=function(e) cat('Current iteration could not be written to ',fileName,'.\n', sep='') )
  } else {
    colnames(tmp) <- c(names(beta),'logLike')
    tryCatch( utils::write.table(tmp,file=fileName, sep=',', row.names=FALSE, append=FALSE),
              error=function(e) cat('Initial iteration could not be written to ',fileName,'.\n', sep='') )
  }
}
