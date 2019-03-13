#' Validates data
#'
#' Checks consistency of the database with \code{apollo_control}, orders it by indivID, and adds an internal ID variable (\code{apollo_sequence})
#'
#' This function should be called after calling apollo_validateControl.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    See \code{?apollo_validateControl} for details.
#' @param silent Boolean. TRUE to keep the function from printing to the console.
#'               Default is FALSE.
#' @return Data.frame. Validated version of database.
apollo_validateData=function(database,apollo_control, silent){

  if(!(apollo_control$indivID %in% names(database))) stop("Column indicated in indivID in apollo_control not found in database. Use valid column name.")

  database <- database[order(database[,apollo_control$indivID]),]

  if(apollo_control$HB==TRUE){
    if("ID" %in% names(database)) warning("Column ID will be used during HB estimation as ID variable.") else warning("You need a column called 'ID' in the database to use HB.")
  }

  database$apollo_sequence <- 1
  tmp <- table(database[,apollo_control$indivID])
  if(length(tmp)>1){
    for(i in 1:length(tmp)) database[database[,apollo_control$indivID]==names(tmp)[i],"apollo_sequence"] <- 1:tmp[i]
  }

  if(!silent) cat("All checks on data completed.\n")
  return(database)

}
