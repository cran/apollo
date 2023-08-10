#' Validates and expands rows if necessary.
#'
#' @param rows Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}. Set to \code{"all"} by default if omitted.
#' @param componentName Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @export
aux_validateRows <- function(rows, componentName=NULL, apollo_inputs=NULL){
  # Error message
  txt1 <- "SYNTAX ISSUE - The \"rows\" argument"
  txt2 <- paste("for model component", componentName)
  txt3 <- paste("needs to be either \"all\" or a vector with as many entries", 
                "as observations in \"database\", where each entry is either", 
                "a logical or 0/1 value.")
  if(!is.null(componentName)){
    txt <- paste(txt1, txt2, txt3)
  } else txt <- paste(txt1, txt3)
  
  # N in database
  test <- is.list(apollo_inputs) || is.environment(apollo_inputs)
  test <- test && !is.null(apollo_inputs$database)
  test <- test && is.data.frame(apollo_inputs$database)
  if(test) N <- nrow(apollo_inputs$database)
  if(!test) stop("INTERNAL ISSUE - Cannot access database")
  
  # Validate rows
  if(is.null(rows)) rows <- "all"
  test <- is.vector(rows)
  if(!test) stop(txt)
  if(is.character(rows)){ # expand if it's "all"
    test <- length(rows)==1 && rows=="all"
    if(test) rows <- return(rep(TRUE, N)) else stop(txt)
  }
  test <- is.logical(rows) || (is.numeric(rows) && all(rows %in% 0:1))
  test <- test && length(rows)==N
  if(!test) stop(txt)
  return(rows)
}