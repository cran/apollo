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
#' @export
apollo_validateData=function(database, apollo_control, silent){
  
  if(!(apollo_control$indivID %in% names(database))) stop("Column indicated in indivID in apollo_control not found in database. Use valid column name!")
  
  if(apollo_control$panelData && length(unique(database[,apollo_control$indivID]))==nrow(database)) stop("Setting panelData in apollo_control is TRUE despite only having one observation per individual in the data!")
  
  # Drop unused levels from factor variables
  database = droplevels.data.frame(database)
  
  # Sort by id
  ##if(apollo_control$panelData) database <- database[order(database[,apollo_control$indivID]),]
  IDs=database[,apollo_control$indivID]
  flagContiguous=TRUE
  for(i in unique(IDs)){
    range=which(IDs==i)
    if(max(range)-min(range)+1!=length(range)) flagContiguous=FALSE
  }
  if(!flagContiguous) stop("All rows for the same individual should be next to each other in the data!")
  
  
  # Rename indivID to "ID" for HB
  if(apollo_control$HB==TRUE){
    if(("ID" %in% names(database))&(apollo_control$indivID!="ID")) apollo_print("Apollo found a column called ID in your database. RSGHB will use this as individual id during HB estimation.") 
    if(!("ID" %in% names(database))) database$ID = database[,apollo_control$indivID]
  }
  
  ### Create scenario id (obs index that resets for each indiv)
  id <- database[,apollo_control$indivID]
  if(!is.numeric(id)) id <- as.numeric(as.factor(id))
  database$apollo_sequence <- sequence(rle(id)$lengths)
  
  # Check existence of weights
  if(!is.null(apollo_control$weights)){
    if(!(apollo_control$weights %in% names(database))) stop("Column ", apollo_control$weights, " not found in database, despite being defined as weights in apollo_control.")
  }
  
  if(!silent) apollo_print("All checks on database completed.")
  return(database)
  
}
