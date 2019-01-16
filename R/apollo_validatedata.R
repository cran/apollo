#' Validates data
#'
#' Checks consistency of the database with \code{apollo_control}, orders it by indivID, and adds an internal ID variable (\code{apollo_ID})
#'
#' This function should be called after calling apollo_validatecontrol.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code.
#'                    See \code{?apollo_validatecontrol} for details.
#' @return Data.frame. Validated version of database.
apollo_validatedata=function(database,apollo_control){

  if(!(apollo_control$indivID %in% names(database))) stop("indivID in apollo_control not found in database. Use valid column name.")

  database <- database[order(database[,apollo_control$indivID]),]

  if(apollo_control$HB==TRUE){
    #if("ID" %in% names(database)) warning("The column called 'ID' was overwritten in the database.")
    #database$ID=database[,apollo_control$indivID]
    if("ID" %in% names(database)) warning("Column ID will be used during HB estimation as individual's id.") else warning("You need a column called 'ID' in the database to use HB.")
  }

  if("apollo_ID" %in% names(database)) warning("The column called 'apollo_ID' was overwritten in the database.")
  database$apollo_ID <- 1
  tmp <- table(database[,apollo_control$indivID])
  if(length(tmp)>1){
    for(i in 1:length(tmp)) database[database[,apollo_control$indivID]==as.numeric(names(tmp)[i]),"apollo_ID"] <- 1:tmp[i]
  }

  cat("All checks on data completed.\n")
  return(database)

}
