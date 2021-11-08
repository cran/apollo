#' Loads model from file
#' 
#' Loads an estimated model object from a file in the current working directory.
#' 
#' This function looks for a file named modelName_model.rds in the working directory, loads the object contained in it, and returns it.
#' @param modelName Character. Name of the model to load.
#' @return A model object.
#' @export
apollo_loadModel <- function(modelName){
  # Validate model name
  if(!(length(modelName)==1 && is.character(modelName))) stop("Argument 'modelName' must be a character variable.")
  
  # Fetch output directory
  outputDirectory <- ''
  tmp  <- tryCatch(get('apollo_control', envir=parent.frame(), inherits=FALSE), error=function(e) FALSE)
  test <- is.list(tmp) && !is.null(tmp$outputDirectory) && is.character(tmp$outputDirectory)
  if(test) outputDirectory <- tmp$outputDirectory else {
    tmp <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE), error=function(e) FALSE)
    test <- is.list(tmp) && !is.null(tmp$apollo_control) && !is.null(tmp$apollo_control$outputDirectory)
    test <- test && is.character(tmp$apollo_control$outputDirectory)
    if(test) outputDirectory <- tmp$apollo_control$outputDirectory
  }
  if(outputDirectory=='') outputDirectory <- getwd()
  test <- !(substr(outputDirectory, nchar(outputDirectory), nchar(outputDirectory)) %in% c('/','\\'))
  if(test) outputDirectory <- paste0(outputDirectory, '/')
  
  # Look for modelname in outputdirectory and working directory, if not found
  fileName <- paste0(outputDirectory, modelName, "_model.rds")
  if(!file.exists(fileName)) fileName <- paste0(getwd(), '/', modelName, "_model.rds")
  model <- tryCatch(readRDS(fileName), error=function(e) stop("Cannot find or open ", fileName) )
  cat("\nSuccessfully loaded ",fileName)
  return(model)
}