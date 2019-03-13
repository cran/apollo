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
#'                     \item{coding}{Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{V}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.}
#'                     \item{rows}{Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded. It must have length equal to the length of argument \code{choiceVar}. Default value is \code{"all"}, meaning all rows are considered in the calculation.}
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
  V=ol_settings[["V"]]
  tau=ol_settings[["tau"]]
  coding=ol_settings[["coding"]]
  rows=ol_settings[["rows"]]

  if(functionality=="validate"){
    apollo_control <- tryCatch(get("apollo_control", parent.frame(), inherits=FALSE),
                            error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))

    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeOrdered))

    if(!apollo_control$noValidation){
      if(!is.vector(tau)) stop("Thresholds for ordered logit needs to be a vector! (no random components allowed).")
      if( !is.null(coding) && (length(tau)+1)!=length(coding) ) stop("Threshold vector length +1 does not match number of elements in argument 'coding'.")
      if(is.null(coding)) coding <- 1:(length(tau)+1)
      if( (length(tau)+1)!=length(unique(outcomeOrdered[rows])) ) stop("Threshold vector length +1 does not match number of unique elements in dependent 'outcomeOrdered'.")
    }

    if(!apollo_control$noDiagnostics){
      choicematrix <- t(as.matrix(table(outcomeOrdered[rows])))
      choicematrix <- rbind(choicematrix, choicematrix[1,]/sum(rows)*100)
      rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")
      print(round(choicematrix,0))
    }

    return(TRUE)
  }

  if(functionality=="zero_LL"){
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeOrdered))

    P <- rep(1/(length(tau)+1),length(outcomeOrdered))
    P[!rows] <- 1
    return(P)
  }

  if(functionality %in% c("estimate","conditionals")){

    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeOrdered))
    if(is.null(coding)) coding <- 1:(length(tau)+1)
    if(functionality=="raw" && length(outcomeOrdered)==1 && is.na(outcomeOrdered)) outcomeOrdered = coding[1]
    outcomeOrdered[!rows]=coding[1]

    map <- stats::setNames(1:length(coding), coding)
    outcomeOrdered2 <- map[as.character(outcomeOrdered)]

    tau <- c(-Inf,tau,Inf)

    p <- 1/(1 + exp(V-tau[outcomeOrdered2+1])) - 1/(1 + exp(V-tau[outcomeOrdered2]))

    if(is.vector(p)) p[!rows]  <- 1
    if(is.matrix(p)) p[!rows,] <- 1
    if(is.array(p) && length(dim(p))==3) p[!rows,,] <- 1

    return(p)
  }

  if(functionality %in% c("prediction", "raw")){

    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeOrdered))
    outcomeOrdered[!rows] <- 1
    tau <- c(-Inf,tau,Inf)

    p = list()
    j=1
    while(j<= length(tau)-1){
      p[[j]] = 1/(1 + exp(V-tau[j+1])) - 1/(1 + exp(V-tau[j]))
      j = j+1
    }

    if(is.null(coding)) coding <- 1:(length(tau)-1)
    names(p) <- coding
	
    if(!(length(outcomeOrdered)==1 && is.na(outcomeOrdered))){
      map <- stats::setNames(1:length(coding), coding)
      outcomeOrdered2  <- map[as.character(outcomeOrdered)]
      p[["chosen"]] <- 1/(1 + exp(V-tau[outcomeOrdered2+1])) - 1/(1 + exp(V-tau[outcomeOrdered2]))
    }
	
    p <- lapply(p, function(x) {
      tmp <- ifelse(functionality=="prediction",NA,1)
      if(is.vector(x)) x[!rows]  <- tmp
      if(is.matrix(x)) x[!rows,] <- tmp
      if(is.array(x) & length(dim(x))==3) x[!rows,,] <- tmp
      return(x)
    })

    return(p)
  }

  if(functionality=="output"){

    p <- apollo_ol(ol_settings, functionality="estimate")

    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(outcomeOrdered))
    choicematrix <- t(as.matrix(table(outcomeOrdered[rows])))
    choicematrix <- rbind(choicematrix, choicematrix[1,]/sum(rows)*100)
    rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")

    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
                                error=function(e){
                                  cat("apollo_ol could not retrieve apollo_control. No diagnostics in output.\n")
                                  return(NA)
                                } )
    if(!(length(apollo_control)==1 && is.na(apollo_control))){
      fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
      fileName <- file.path(tempdir(),fileName)
      fileConn <- tryCatch( file(fileName, open="at"),
                            error=function(e){
                              cat('apollo_ol could not write diagnostics to temporary file. No diagnostics in output.\n')
                              return(NA)
                            })
      if(!anyNA(fileConn)){
        sink(fileConn)
        on.exit({if(sink.number()>0) sink(); close(fileConn)})
        if(apollo_control$noDiagnostics==FALSE){
          cat('Overview of choices for OL model component:\n')
          print(round(choicematrix,0))
          cat('\n')}
      }
    }

    return(p)
  }


}
