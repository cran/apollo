#' Calculates multinomial logit probabilities
#'
#' Calculates probabilities of a multinomial logit model.
#'
#' @param mnl_settings List of inputs of the MNL model. It should contain the following.
#'                     \itemize{
#'                       \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item choiceVar: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item V: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                     }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate": Used for model estimation.
#'                        \item "prediction": Used for model predictions.
#'                        \item "validate": Used for validating input.
#'                        \item "zero_LL": Used for calculating null likelihood.
#'                        \item "conditionals": Used for calculating conditionals.
#'                        \item "output": Used for preparing output after model estimation.
#'                        \item "raw": Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @examples
#' ### Load data
#' data(apollo_modeChoiceData)
#' database <- apollo_modeChoiceData
#' rm(apollo_modeChoiceData)
#'
#' ### Parameters
#' b = list(asc_1=0, asc_2=0, asc_3=0, asc_4=0, tt=0, tc=0, acc=0)
#'
#' ### List of utilities
#' V = list()
#' V[['car' ]] = b$asc_1 + b$tt*database$time_car  + b$tc*database$cost_car
#' V[['bus' ]] = b$asc_2 + b$tt*database$time_bus  + b$tc*database$cost_bus  +
#'               b$acc*database$access_bus
#' V[['air' ]] = b$asc_3 + b$tt*database$time_air  + b$tc*database$cost_air  +
#'               b$acc*database$access_air
#' V[['rail']] = b$asc_4 + b$tt*database$time_rail + b$tc*database$cost_rail +
#'               b$acc*database$access_rail
#'
#' ### MNL settings
#' mnl_settings <- list(
#'    alternatives = c(car=1, bus=2, air=3, rail=4),
#'    avail        = list(car=database$av_car, bus=database$av_bus,
#'                        air=database$av_air, rail=database$av_rail),
#'    choiceVar    = database$choice,
#'    V            = V
#' )
#'
#' ### Compute choice probabilities using MNL model
#' apollo_mnl(mnl_settings, functionality="estimate")
#' @export
apollo_mnl <- function(mnl_settings, functionality){
  if(is.null(mnl_settings[["alternatives"]])) stop("The mnl_settings list needs to include an object called \"alternatives\"!")
  if(is.null(mnl_settings[["avail"]])) stop("The mnl_settings list needs to include an object called \"avail\"!")
  if(is.null(mnl_settings[["choiceVar"]])) stop("The mnl_settings list needs to include an object called \"choiceVar\"!")
  if(is.null(mnl_settings[["V"]])) stop("The mnl_settings list needs to include an object called \"V\"!")
  if(is.null(mnl_settings[["rows"]])) mnl_settings[["rows"]]="all"

  alternatives = mnl_settings[["alternatives"]]
  avail        = mnl_settings[["avail"]]
  choiceVar    = mnl_settings[["choiceVar"]]
  V            = mnl_settings[["V"]]
  rows         = mnl_settings[["rows"]]

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    # Store useful values
    apollo_control <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_control,
                            error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # check rows statement
    if(length(rows)!=length(choiceVar)) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of the rows in the data!")
    
    # Check that alternatives are named in altcodes and V
    if(is.null(altnames) || is.null(altcodes) || is.null(names(V))) stop("Alternatives must be named, both in 'altnernatives' and 'V'.")
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
      #warning("Full availability of alternatives assumed for the MNL component.")
    }
    # Reorder V to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]

    if(apollo_control$noValidation==FALSE){
      # Check that there are at least two alternatives
      if(nAlts<2) stop("MNL requires at least two alternatives")

      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")
      
      choiceLabs <- unique(choiceVar[rows])
      if(!all(altnames %in% names(V))) stop("Alternatives names in \"alternatives\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternatives names in \"alternatives\" do not match those in \"avail\".")

      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop('A value in "choiceVar" column is not included in "alternatives".')

      # check that nothing unavailable is chosen
      chosenunavail=0
      j=1
      while(j <= length(altnames)){
        if(sum((choiceVar==altcodes[j])*(avail[[j]]==0)*rows)) chosenunavail=1
        j=j+1
      }
      if(chosenunavail==1) stop("Some alternative(s) chosen despite being listed as unavailable\n")

      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")

      #cat("\nAll checks passed for MNL model component\n")

    }

    if(apollo_control$noDiagnostics==FALSE){
      

      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)

      # Construct summary table of availabilities and market share
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      #choicematrix[1,] = colSums(rows*matrix(unlist(avail), ncol = length(avail)))
      choicematrix[1,] = unlist(lapply(avail, function(x) sum(x[rows])))
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j] = sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = altnames
      
      content <- list(round(choicematrix,2),
                      ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
                      ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""))
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided to 'apollo_mnl' (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for MNL model component:", content, apolloLog)
      
      #cat('Overview of choices for MNL model component:\n')
      #print(round(choicematrix,2))
      #cat("\n")
      #if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      #if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
    }

    return(invisible(TRUE))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]

    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }

    if(!anyNA(avail)) if(all(altnames != names(avail))) avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(matrix(unlist(avail), ncol = length(avail))) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    P[!rows] = 1
    return(P)
  }

  # ############################################################ #
  #### functionality="estimate/prediction/conditionals/raw" ####
  # ############################################################ #

  if(functionality %in% c("estimate","prediction","conditionals","raw")){

    # Fix choiceVar if "raw" and choiceVar==NA
    choiceNA = FALSE
    if(length(choiceVar)==1 && is.na(choiceVar)){
      choiceVar = alternatives[1]
      choiceNA = TRUE
    }

    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]

    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }

    # Set utility of unavailable alternatives to 0 (to avoid numerical issues)
    for(a in 1:nAlts){
      isSca <- length(V[[a]])==1
      isVec <- is.vector(V[[a]]) & !isSca
      isMat <- is.matrix(V[[a]])
      isCub <- is.array(V[[a]]) && !isMat && length(dim(V[[a]]))==3
      if(isVec) V[[a]][!avail[[a]]]   <- 0
      if(isMat) V[[a]][!avail[[a]],]  <- 0
      if(isCub) V[[a]][!avail[[a]],,] <- 0
    }

    ### Reorder V and avail to match altnames order, if necessary
    if(all(altnames != names(V))) V <- V[altnames]
    if(all(altnames != names(avail))) avail <- avail[altnames]

    # extract chosen utility
    chosenV <- Reduce('+',
                      lapply(as.list(1:nAlts),
                             FUN=function(i) (choiceVar==altcodes[i])*V[[altnames[i]]])
    )

    # substract chosen utility from all others for numerical stability
    # uses lapply to loop over individual matrices in list of utilities
    V = lapply(X=V, FUN=function(v) v-chosenV)

    # consider availabilities once before exponentiating (avoids issues if unavailable alternatives have attributes at 999)
    V <- mapply('*', V, avail, SIMPLIFY = FALSE)

    # exponentiate utilities
    eV = lapply(X=V, FUN=exp)

    # consider availabilities (it assumes eV and avail are in the same order)
    eV <- mapply('*', eV, avail, SIMPLIFY = FALSE)

    # calculate the denominator of the logit probability expression
    # now simply the addition of all exponentiated utilities
    denom = Reduce('+',eV)

    # calculate probabilities
    # if probabilities for all alternatives are requested, then P is a list
    # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
    if((functionality=="prediction")|(functionality=="raw")){
      P <- lapply(eV, function(ev) ev/denom)
      if(!choiceNA) P[["chosen"]]= 1/denom
      P <- lapply(P, function(p) {
        # change likelihood of excluded columns
        if(is.vector(p)) p[!rows]  <- ifelse(functionality=="prediction",NA,1)
        if(is.matrix(p)) p[!rows,] <- ifelse(functionality=="prediction",NA,1)
        if(is.array(p) & length(dim(p))==3) p[!rows,,] <- ifelse(functionality=="prediction",NA,1)
        return(p)
      })
    } else {
      # change likelihood of excluded columns to 1
      P <- 1/denom
      if(is.vector(P)) P[!rows]  <- 1
      if(is.matrix(P)) P[!rows,] <- 1
      if(is.array(P) && length(dim(P))==3) P[!rows,,] <- 1
    }

    return(P)
  }

  # ############################## #
  #### functionality="output" ####
  # ############################## #

  if(functionality=="output"){

    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]

    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }

    ### Reorder V and avail to match altnames order, if necessary
    if(all(altnames != names(V))) V <- V[altnames]
    if(all(altnames != names(avail))) avail <- avail[altnames]

    # extract chosen utility
    chosenV <- Reduce('+',
                      lapply(as.list(1:nAlts),
                             FUN=function(i) (choiceVar==altcodes[i])*V[[altnames[i]]])
    )

    # substract chosen utility from all others for numerical stability
    # uses lapply to loop over individual matrices in list of utilities
    V = lapply(X=V, FUN=function(v) v-chosenV)

    # consider availabilities once before exponentiating (avoids issues if unavailable alternatives have attributes at 999)
    V <- mapply('*', V, avail, SIMPLIFY = FALSE)

    # exponentiate utilities
    eV = lapply(X=V, FUN=exp)

    # consider availabilities (it assumes eV and avail are in the same order)
    eV <- mapply('*', eV, avail, SIMPLIFY = FALSE)

    # calculate the denominator of the logit probability expression
    # now simply the addition of all exponentiated utilities
    denom = Reduce('+',eV)

    # calculate probabilities
    P <- 1/denom

    # change likelihood of excluded columns to 1
    if(is.vector(P)) P[!rows]  <- 1
    if(is.matrix(P)) P[!rows,] <- 1
    if(is.array(P) && length(dim(P))==3) P[!rows,,] <- 1

    # turn scalar availabilities into vectors
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)

    # Construct summary table of availabilities and market share
    choicematrix = matrix(0,nrow=4,ncol=length(altnames))
    choicematrix[1,] = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
    choicematrix[1,] = unlist(lapply(avail, function(x) sum(x[rows])))
    j=1
    while(j<= length(altnames)){
      choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
      j=j+1
    }
    choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
    choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
    choicematrix[4,!is.finite(choicematrix[4,])] <- 0
    rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
    colnames(choicematrix) = altnames
    
    content <- list(round(choicematrix,2),
                    ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
                    ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""))
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided to 'apollo_mnl' (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for MNL model component:", content, apolloLog)
    
    ## write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
    #apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
    #                            error=function(e){
    #                              cat("apollo_mnl could not retrieve apollo_control. No diagnostics in output.\n")
    #                              return(NA)
    #                            } )
    #if(!(length(apollo_control)==1 && is.na(apollo_control))){
    #  fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    #  fileName <- file.path(tempdir(),fileName)
    #  fileConn <- tryCatch( file(fileName, open="at"),
    #                        error=function(e){
    #                          cat('apollo_mnl could not write diagnostics to temporary file. No diagnostics in output.\n')
    #                          return(NA)
    #                        })
    #  if(!anyNA(fileConn)){
    #    sink(fileConn)
    #    on.exit({if(sink.number()>0) sink(); close(fileConn)})
    #    if(apollo_control$noDiagnostics==FALSE){
    #      cat('Overview of choices for MNL model component:\n')
    #      print(round(choicematrix,2))
    #      cat("\n")}
    #    if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
    #    if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
    #  }
    #}


    return(P)
  }
}