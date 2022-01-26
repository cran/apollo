#' Writes the vector [beta,ll] to a file called \code{modelname_iterations.csv}
#' @param beta vector of parameters to be written (including fixed ones).
#' @param ll scalar representing the log-likelihood of the whole model.
#' @param modelName Character. Name of the model.
#' @return Nothing.
#' @export
apollo_writeTheta <- function(beta, ll, modelName){
  ### Fetch apollo_inputs
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(1), inherits=TRUE), error=function(e) NULL)
  
  ### Scale parameters
  test <- !is.null(apollo_inputs) && is.list(apollo_inputs) && !is.null(apollo_inputs$apollo_scaling)
  test <- test && is.numeric(apollo_inputs$apollo_scaling) && is.vector(apollo_inputs$apollo_scaling)
  test <- test && !is.null(names(apollo_inputs$apollo_scaling))
  if(test) beta[names(apollo_inputs$apollo_scaling)] <- beta[names(apollo_inputs$apollo_scaling)]*apollo_inputs$apollo_scaling
  
  ### Set output directory
  test <- !is.null(apollo_inputs) && is.list(apollo_inputs) && !is.null(apollo_inputs$apollo_control)
  test <- test && is.list(apollo_inputs$apollo_control) && !is.null(apollo_inputs$apollo_control$outputDirectory)
  test <- test && is.character(apollo_inputs$apollo_control$outputDirectory)
  test <- test && length(apollo_inputs$apollo_control$outputDirectory)==1
  if(test) outputDirectory <- apollo_inputs$apollo_control$outputDirectory else outputDirectory=paste0(getwd(),'/')
  
  ### Sort parameters in the same order than apollo_beta
  test <- !is.null(apollo_inputs) && is.list(apollo_inputs)
  test <- test && !is.null(apollo_inputs$apollo_beta_names) && is.vector(apollo_inputs$apollo_beta_names)
  test <- test && is.character(apollo_inputs$apollo_beta_names) 
  test <- test && all(apollo_inputs$apollo_beta_names %in% names(beta))
  if(test) beta <- beta[apollo_inputs$apollo_beta_names]
  
  ### Initialise
  fileName <- paste0(outputDirectory,modelName,"_iterations.csv")
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
