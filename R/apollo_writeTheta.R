#' Writes the vector [beta,ll] to a file called modelname_iterations.csv
#' @param beta vector of parameters to be written.
#' @param ll scalar representing the loglikelihood of the whole model.
#' @param modelName Character. Name of the model.
#' @return Nothing.
#' @export
apollo_writeTheta <- function(beta, ll, modelName){
  
  ### Fetch apollo_control
  # Try geting it from apollo_inputs
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(1), inherits=FALSE),
                            error=function(e) NULL)
  test <- !is.null(apollo_inputs) && !is.null(apollo_inputs$apollo_control)
  if(test) apollo_control <- apollo_inputs$apollo_control else apollo_control <- NULL
  rm(apollo_inputs)
  # If apollo_inputs failed, try getting it from the calling environment
  if(is.null(apollo_control)) apollo_control <- tryCatch(get('apollo_control', envir=parent.frame(1), inherits=FALSE),
                                                         error=function(e) NULL)
  # If could not be fetched, set to working directory
  apollo_control <- list(outputDirectory=paste0(getwd(),'/'))
  
  
  ### Initialise
  fileName <- paste0(apollo_inputs$apollo_control$outputDirectory,modelName,"_iterations.csv")
  tmp <- matrix(c(beta,ll ),nrow=1)
  
  if(file.exists(fileName)){
    # If file already exists, append
    tryCatch( utils::write.table(tmp,file=fileName, append=TRUE, sep=',', col.names=FALSE, row.names=FALSE),
              error=function(e) cat('Current iteration could not be written to ',fileName,'.\n', sep='') )
  } else {
    # If file does not exist, write afresh with headings
    colnames(tmp) <- c(names(beta),'logLike')
    tryCatch( utils::write.table(tmp,file=fileName, sep=',', row.names=FALSE, append=FALSE),
              error=function(e) cat('Initial iteration could not be written to ',fileName,'.\n', sep='') )
  }
}
