#' Validates data
#'
#' Checks consistency of the database with \code{apollo_control}, sorts it by indivID, and adds an internal ID variable (\code{apollo_sequence})
#'
#' This function should be called after calling apollo_validateControl.
#' Observations are sorted only if apollo_control$panelData=TRUE.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    See \code{?apollo_validateControl} for details.
#' @param silent Boolean. TRUE to keep the function from printing to the console.
#'               Default is FALSE.
#' @return Data.frame. Validated version of database.
apollo_validateData=function(database, apollo_control, silent){
  
  if(!(apollo_control$indivID %in% names(database))) stop("Column indicated in indivID in apollo_control not found in database. Use valid column name.")
  
  # Drop unused levels from factor variables
  database = droplevels.data.frame(database)
  
  # Sort by id
  if(apollo_control$panelData) database <- database[order(database[,apollo_control$indivID]),]
  
  # Rename indivID to "ID" for HB
  if(apollo_control$HB==TRUE){
    if(("ID" %in% names(database))&(apollo_control$indivID!="ID")) cat("\nColumn ID will be used during HB estimation as ID variable.") 
    if(!("ID" %in% names(database))) database$ID = database[,apollo_control$indivID]
  }
  
  ### Create scenario id (obs index that resets for each indiv)
  id <- database[,apollo_control$indivID]
  if(!is.numeric(id)) id <- as.numeric(as.factor(id))
  database$apollo_sequence <- sequence(rle(id)$lengths)
  
  # Check existance of weights
  if(!is.null(apollo_control$weights)){
    if(!(apollo_control$weights %in% names(database))) stop("Column ", apollo_control$weights, " not found in database, despite being defined as weights in 'apollo_control'.")
  }
  
  if(!silent) cat("All checks on data completed.\n")
  return(database)
  
}
