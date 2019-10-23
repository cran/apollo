#' Calculates the probability of an ordered logit model
#'
#' Calculates the probabilities of an ordered logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' This function estimates an ordered logit model of the type:
#' y* = V + epsilon
#' outcomeOrdered =  1 if   -Inf < y* < tau[1]
#'      2 if tau[1] < y* < tau[2]
#'      ...
#'      maxLvl if tau[length(tau)] < y* < +Inf
#' Where epsilon is distributed standard logistic, and the values 1, 2, ..., maxLvl can be
#' replaces by coding[1], coding[2], ..., coding[maxLvl].
#' The behaviour of the function changes depending on the value of the \code{functionality} argument.
#' @param ol_settings List of settings for the OL model. It should include the following.
#'                   \describe{
#'                     \item{outcomeOrdered}{Numeric vector. Dependant variable. The coding of this variable is assumed to be from 1 to the maximum number of different levels. For example, if the ordered response has three possible values: "never", "sometimes" and "always", then it is assumed that outcomeOrdered contains "1" for "never", "2" for "sometimes", and 3 for "always". If another coding is used, then it should be specified using the \code{coding} argument.}
#'                     \item{V}{Numeric vector. A single explanatory variable (usually a latent variable). Must have the same number of rows as outcomeOrdered.}
#'                     \item{tau}{Numeric vector. Thresholds. As many as number of different levels in the dependent variable - 1. Extreme thresholds are fixed at -inf and +inf. No mixing allowed in thresholds.}
#'                     \item{coding}{Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{outcomeOrdered}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.}
#'                     \item{rows}{Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded. It must have length equal to the length of argument \code{outcomeOrdered}. Default value is \code{"all"}, meaning all rows are considered in the calculation.}
#'                   }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \describe{
#'                        \item{"estimate"}{Used for model estimation.}
#'                        \item{"prediction"}{Used for model predictions.}
#'                        \item{"validate"}{Used for validating input.}
#'                        \item{"zero_LL"}{Used for calculating null likelihood.}
#'                        \item{"conditionals"}{Used for calculating conditionals.}
#'                        \item{"output"}{Used for preparing output after model estimation.}
#'                        \item{"raw"}{Used for debugging.}
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \describe{
#'           \item{"estimate"}{vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.}
#'           \item{"prediction"}{List. Returns a list with the probabilities for each possible outcome.}
#'           \item{"validate"}{Boolean. Returns TRUE if all tests are passed.}
#'           \item{"zero_LL"}{vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.}
#'           \item{"conditionals"}{Same as "prediction".}
#'           \item{"output"}{Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).}
#'           \item{"raw"}{Same as "prediction".}
#'         }
#' @export
apollo_ol  <- function(ol_settings, functionality){
  if(is.null(ol_settings[["outcomeOrdered"]])) stop("The ol_settings list needs to include an object called \"outcomeOrdered\"!")
  if(is.null(ol_settings[["V"]])) stop("The ol_settings list needs to include an object called \"V\"!")
  if(is.null(ol_settings[["tau"]])) stop("The ol_settings list needs to include an object called \"tau\"!")
  if(is.null(ol_settings[["coding"]])) ol_settings[["coding"]]=NULL
  if(is.null(ol_settings[["rows"]])) ol_settings[["rows"]]="all"

  outcomeOrdered=ol_settings[["outcomeOrdered"]]
  V      = ol_settings[["V"]]
  tau    = ol_settings[["tau"]]
  coding = ol_settings[["coding"]]
  rows   = ol_settings[["rows"]]
  nObs <- tryCatch(nrow( get("apollo_inputs", parent.frame(), inherits=FALSE)$database ),
                   error=function(e){
                     lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenC <- length(outcomeOrdered)
                     return(max(lenV, lenC))
                   })
  
  ### Format checks
  # outcomeOrdered
  test <- is.vector(outcomeOrdered) && (length(outcomeOrdered)==nObs || length(outcomeOrdered)==1)
  if(!test) stop("The \"outcomeOrdered\" argument needs to be a scalar or a vector with one entry per observation in the \"database\"")
  # V
  test <- is.numeric(V) && ( is.vector(V) || (is.array(V) && length(dim(V))<=3) )
  if(!test) stop("Argument \"V\" must be a numeric vector, matrix or 3-dim array.")
  test <- (is.vector(V) && length(V) %in% c(1,nObs)) || (is.array(V) && dim(V)[1] %in% c(1,nObs))
  if(!test) stop("Argument \"V\" must be have 1 row, or as many as the \"database\".")
  # rows
  test <- is.vector(rows) && ( (is.logical(rows) && length(rows)==nObs) || (length(rows)==1 && rows=="all") )
  if(!test) stop("The \"rows\" argument needs to be \"all\" or a vector of boolean statements with one entry per observation in the \"database\"")
  # functionality
  test <- functionality %in% c("estimate","prediction","validate","zero_LL","conditionals","output","raw")
  if(!test) stop("Non-permissable setting for \"functionality\"")
  # coding
  test <- is.null(coding) || (is.vector(coding))
  if(!test) stop("Argument 'coding', if provided, must be a vector.")
  
  ### Get number of observations
  nObs <- tryCatch(nrow( get("apollo_inputs", parent.frame(), inherits=TRUE)$database ),
                   error=function(e){
                     lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenC <- length(outcomeOrdered)
                     return(max(lenV, lenC))
                   })
  
  ### Validate input
  #if(!is.vector(outcomeOrdered)) stop("Argument 'outcomeOrdered' must be a vector of length nObs.")
  #if(!is.numeric(V) || !(is.array(V) | is.vector(V))) stop("Argument 'V' must be a numeric vector or array.")
  #if(!is.null(coding) && (!is.vector(coding) | is.null(names(coding)))) stop("Argument 'coding', if provided, must be a named vector.")
  ### CORRECTED NEXT LINE
  ### if(rows!="all" & ( length(rows)!=nObs | !is.logical(rows) )) stop("Argument 'rows', if provided, must be \"all\" or a vector of boolean statements of length nObs.")
  #if(!( length(rows)==1 && (rows=="all")) & ( length(rows)!=nObs | !is.logical(rows) )) stop("Argument 'rows', if provided, must be \"all\" or a vector of boolean statements of length nObs.")
  
  ### Filter rows
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, length(nObs))
  if(any(!rows)){
    outcomeOrdered <- outcomeOrdered[rows]
    V <- apollo_keepRows(V, rows)
  }
  nObs <- max(ifelse(is.array(V), dim(V)[1], length(V)) , length(outcomeOrdered))
  

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    apollo_control <- tryCatch(get("apollo_control", parent.frame(), inherits=FALSE),
                            error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    
    if(!apollo_control$noValidation){
      if(!is.vector(tau)) stop("Thresholds for ordered logit needs to be a vector (no random components allowed)!")
      values_present=unique(outcomeOrdered)
      if(is.null(coding)){
        coding <- 1:(length(tau)+1)  
        cat("\nNo coding provided for ordered logit, so assuming\n outcomeOrdered goes from 1 to",max(coding),"\n")
      }
      if(!(all(values_present %in% coding ))) stop("Some levels specified in outcomeOrdered do not exist in coding!")
      if(!(all(coding %in% values_present ))) stop("Some levels specified in coding do not exist in outcomeOrdered!")
      if( (length(tau)+1)!=length(coding) ) stop("Threshold vector length +1 does not match number of elements in argument 'coding'.")
    }
    
    if(!apollo_control$noDiagnostics){
      choicematrix <- t(as.matrix(table(outcomeOrdered)))
      choicematrix <- rbind(choicematrix, choicematrix[1,]/nObs*100)
      rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for OL model component:", round(choicematrix,0), apolloLog)
    }

    return(TRUE)
  }

  # ########################################################## #
  #### functionality="zero_LL"                              ####
  # ########################################################## #

  if(functionality=="zero_LL"){
    P <- rep(NA, nObs)
    if(any(!rows)) P <- apollo_insertRows(P, rows, 1)
    return(P)
  }

  # ########################################################## #
  #### functionality="estimate | conditionals | raw"        ####
  # ########################################################## #

  if(functionality %in% c("estimate","conditionals")){

    if(is.null(coding)) coding <- 1:(length(tau)+1)
    map <- stats::setNames(1:length(coding), coding)
    outcomeOrdered2 <- map[as.character(outcomeOrdered)]

    tau <- c(-Inf,tau,Inf)

    p <- 1/(1 + exp(V-tau[outcomeOrdered2+1])) - 1/(1 + exp(V-tau[outcomeOrdered2]))

    if(any(!rows)) p <- apollo_insertRows(p, rows, 1) # insert excluded rows with value 1
    
    return(p)
  }

  # ########################################################## #
  #### functionality="prediction"                           ####
  # ########################################################## #
  if(functionality %in% c("prediction", "raw")){
    
    if(is.null(coding)) coding <- 1:(length(tau)+1)
    tau <- c(-Inf,tau,Inf)
    
    p = list()
    for(j in 1:(length(tau)-1)) p[[j]] = 1/(1 + exp(V-tau[j+1])) - 1/(1 + exp(V-tau[j]))
    
    if(is.null(coding)) coding <- 1:(length(tau)-1)
    names(p) <- coding
    
    # insert excluded rows with value NA
	  if(any(!rows)) p <- lapply(p, apollo_insertRows, r=rows, val=NA)
    
    # Add chosen alternative (unless outcomeOrdered is NA)
    if(!(length(outcomeOrdered)==1 && is.na(outcomeOrdered))){
      p[["chosen"]] <- apollo_ol(ol_settings, functionality="estimate")
    }

    return(p)
  }

  # ############################# #
  #### functionality="output"  ####
  # ############################# #
  if(functionality=="output"){

    p <- apollo_ol(ol_settings, functionality="estimate")
    
    choicematrix <- t(as.matrix(table(outcomeOrdered)))
    choicematrix <- rbind(choicematrix, choicematrix[1,]/nObs*100)
    rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")
    
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for OL model component:", round(choicematrix,0), apolloLog)

    return(p)
  }


}
