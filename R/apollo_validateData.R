#' Validates data
#'
#' Checks consistency of the database with \code{apollo_control}, sorts it by 
#' indivID, and adds an internal ID variable (\code{apollo_sequence})
#'
#' This function should be called after calling \link{apollo_validateControl}.
#' Observations are sorted only if \code{apollo_control$panelData=TRUE}.
#' @param database data.frame. Data used by model.
#' @param apollo_control List. Options controlling the running of the code. See 
#'                       \link{apollo_validateInputs}.
#' @param silent Boolean. TRUE to prevent the function from printing to the 
#'               console. Default is FALSE.
#' @return Data.frame. Validated version of database.
#' @export
apollo_validateData=function(database, apollo_control, silent){
  
  if(nrow(database)==0) stop("INPUT ISSUE - database is empty!")
  
  if(any(is.na(database))){
    txt <- paste0("Your database contains some entries that are NA. ", 
                  "This may well be intentional, but be advised that if ", 
                  "these entries are used in your model, the behaviour may ", 
                  "be unexpected.")
    if(!silent) apollo_print(txt, pause=0, type="w")}
  
  test <- apollo_control$indivID %in% names(database)
  if(!test) stop("INPUT ISSUE - Column indicated in indivID in apollo_control not found in ", 
                 "database. Use valid column name!")
  
  test <- apollo_control$panelData
  test <- test && length(unique(database[,apollo_control$indivID]))==nrow(database)
  if(test) stop("INCORRECT FUNCTION/SETTING USE - Setting panelData in apollo_control is TRUE despite only ", 
                "having one observation per individual in the data!")
  
  # Drop unused levels from factor variables
  database = droplevels.data.frame(database)
  
  # Sort by id
  ##if(apollo_control$panelData) database <- database[order(database[,apollo_control$indivID]),]
  IDs = database[,apollo_control$indivID]
  flagContiguous = TRUE
  for(i in unique(IDs)){
    range = which(IDs==i)
    if(max(range)-min(range)+1 != length(range)) flagContiguous=FALSE
  }
  if(!flagContiguous) stop("INPUT ISSUE - All rows for the same individual should be next ", 
                           "to each other in the data!")
  
  
  # Rename indivID to "ID" for HB
  if(apollo_control$HB==TRUE){
    test <- ("ID" %in% names(database)) && apollo_control$indivID!="ID"
    if(test) if(!silent) apollo_print(paste0("Apollo found a column called ID in your ", 
                                             "database. RSGHB will use this as individual ",
                                             "id during HB estimation."), type="i")
    test <- "ID" %in% names(database)
    if(!test) database$ID = database[,apollo_control$indivID]
  }
  
  ### Create scenario id (obs index that resets for each indiv)
  id <- database[,apollo_control$indivID]
  if(!is.numeric(id)) id <- as.numeric(as.factor(id))
  database$apollo_sequence <- sequence(rle(id)$lengths)
  
  # Check existence of weights
  if(!is.null(apollo_control$weights)){
    test <- apollo_control$weights %in% names(database)
    if(!test) stop("INPUT ISSUE - Column ", apollo_control$weights, " not found in the ", 
                   "database despite being defined as weights in ", 
                   "apollo_control.")
  }
  
  # Warn of the presence of "factor" variables
  isFactor <- sapply(database, is.factor)
  if(any(isFactor)){
    #data("apollo_modeChoiceData")
    #database <- apollo_modeChoiceData
    #database$choice <- as.factor(database$choice)
    #database$female <- as.factor(database$female)
    #isFactor <- sapply(database, is.factor)
    txt <- paste0(names(database)[isFactor], collapse="\", \"")
    txt <- paste0("WARNING: Your database contains variable(s) \"", txt, "\" ", 
                  "codified as factors. Apollo does not support factors, ", 
                  "and using them inside apollo_probabilities may lead to ", 
                  "NA values in the loglikelihood. If you want to use ", 
                  "these variables, we recommend manually transforming them ", 
                  "into numeric variables.")
    if(!silent) apollo_print(txt)
  }
  
  
  if(!silent) apollo_print("All checks on database completed.")
  return(database)
  
}
