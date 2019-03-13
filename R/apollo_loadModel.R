#' Loads model from file
#' 
#' Loads an estimated model object from a file in the current working directory.
#' 
#' This function looks for a file named modelName_model.rds in the working directory, loads the object contained in it, and returns it.
#' @param modelName Character. Name of the model to load.
#' @return A model object.
#' @export
apollo_loadModel <- function(modelName){
  if(!(length(modelName)==1 && is.character(modelName))) stop("Argument 'modelName' must be a character variable.")
  
  fileName <- paste0(modelName, "_model.rds")
  model <- tryCatch(readRDS(fileName), error=function(e) stop("Cannot find or open ", fileName) )
  cat("\nSuccessfully loaded ",fileName)
  return(model)
}