#' Calculates multinomial logit probabilities
#'
#' Calculates the probabilities of a multinomial logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' The behaviour of the function depends on the value of the \code{functionality} argument.
#' Calculated probabilities will take the same form than V. The shape of V varies automatically depending
#' on whether draws are used or not, and if used, what kind of draws they are. In general, this should not
#' be a concern for the user, as all re-shapes and dimension matching are handled internally. Only when
#' the user does not use pre-coded likelihood functions (such as \link{apollo_mnl}), and instead programs
#' his/her own likelihood, should the user be concern about the dimensionality of V. For those cases,
#' consider the following guide.
#' \describe{
#'   \item{No mixing}{V should be a vector with as many elements as observations in the database (nObs).}
#'   \item{inter-individual draws only}{V should be a matrix with as many rows as observations (nObs) and as many columns as inter-individual draws (nInterDraws).}
#'   \item{intra-individual draws}{V should be a 3-dimensional array with dimensions nObs x nInterDraws x nIntraDraws. If nInterDraws=0, then dimensions are nObs x 1 x nIntraDraws.}
#' }
#' @param alternatives Named numeric vector. This vector contains the name of the alternatives and the
#'                     value representing them in the \code{choice} argument. All different elements in argument
#'                     \code{choice} must be included in argument \code{alternatives}.
#' @param avail Named list. Availabilities of alternatives, one element per alternative.
#'              Names of elements must match those in argument \code{alternatives}.
#'              Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation.
#'              If all alternatives then the user can just set \code{avail=1}.
#' @param choice Numeric vector. Contains choices for all observations. It will usually be a column from the database.
#' @param V Named list. Utilities, names of elements must match those in argument \code{alternatives}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \describe{
#'                        \item{"estimate"}{Used for model estimation.}
#'                        \item{"prediction"}{Used for model predictions.}
#'                        \item{"validate"}{Used for validating input.}
#'                        \item{"zero_LL"}{Used for getting null likelihood.}
#'                        \item{"conditionals"}{Used for getting conditionals.}
#'                        \item{"output"}{Used for preparing output after model estimation.}
#'                        \item{"raw"}{Used for debugging.}
#'                      }
#' @param rows Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded.
#'             It must have length equal to the length of argument \code{choice}.
#'             Default value is \code{"all"}, meaning all rows are considered in the calculation.
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#' \describe{
#'   \item{"estimate"}{vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.}
#'   \item{"prediction"}{List. Returns a list with the probabilities for all alternatives and additionaly the probability of the chosen alternative.}
#'   \item{"validate"}{Boolean. Returns TRUE if all tests are passed.}
#'   \item{"zero_LL"}{vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.}
#'   \item{"conditionals"}{Same as "estimation".}
#'   \item{"output"}{Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modeloutput}).}
#'   \item{"raw"}{Same as "prediction".}
#' }
#'
#' @examples
#' data(apollo_modeChoiceData)
#' x            = database
#' b            = list(asc_1=0, asc_2=0, asc_3=0, asc_4=0, tt=0, tc=0, acc=0)
#' alternatives = c(car=1, bus=2, air=3, rail=4)
#' avail        = list(car=x$av_car, bus=x$av_bus, air=x$av_air, rail=x$av_rail)
#' choiceVar    = x$choice
#'
#' ### List of utilities
#' V = list()
#' V[['car' ]] = b$asc_1 + b$tt*x$time_car  + b$tc*x$cost_car
#' V[['bus' ]] = b$asc_2 + b$tt*x$time_bus  + b$tc*x$cost_bus  + b$acc*x$access_bus
#' V[['air' ]] = b$asc_3 + b$tt*x$time_air  + b$tc*x$cost_air  + b$acc*x$access_air
#' V[['rail']] = b$asc_4 + b$tt*x$time_rail + b$tc*x$cost_rail + b$acc*x$access_rail
#'
#' ### Compute choice probabilities using MNL model
#' apollo_mnl(alternatives, avail, choiceVar, V)
#'
#' @export
apollo_mnl <- function(alternatives, avail, choice, V, functionality="estimate", rows="all"){

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    # Store useful values
    apollo_control <- tryCatch(get("apollo_control", parent.frame(), inherits=FALSE),
                            error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    nObs  <- length(choice)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choice))
    choice[!rows]=alternatives[1]

    # check rows statement
    if(length(rows)!=length(choice)) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of the rows in the data!")

    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
      warning("Full availability of alternatives assumed for the MNL component.")
    }
    # Reorder V to match altnames order, if necessary
    if(all(altnames != names(V))) V <- V[altnames]
    if(all(altnames != names(avail))) avail <- avail[altnames]

    if(apollo_control$noValidation==FALSE){
      # Check that there are at least two alternatives
      if(nAlts<2) stop("MNL requires at least two alternatives")

      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")

      # Check that labels in choice match those in the utilities and availabilities
      choiceLabs <- unique(choice)
      if(!all(altnames %in% names(V))) stop("Alternative labels in \"altnames\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternative labels in \"altnames\" do not match those in \"avail\".")

      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop("Value in choice column that is not included in altcodes.")

      # check that nothing unavailable is chosen
      chosenunavail=0
      j=1
      while(j <= length(altnames)){
        if(sum((choice==altcodes[j])*(avail[[j]]==0)*rows)) chosenunavail=1
        j=j+1
      }
      if(chosenunavail==1) stop("Some alternative(s) chosen despite being listed as unavailable\n")

      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")

      cat("\nAll checks passed for MNL model component\n")

    }

    if(apollo_control$noDiagnostics==FALSE){
      if(avail_set==TRUE) warning("Availability not provided to 'apollo_mnl' (or some elements are NA). Full availability assumed.")

      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)

      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j]=sum(choice==altcodes[j]) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      rownames(choicematrix) = c("Times available","Times chosen","Percentage of choice overall","Percentage of choice when available")
      colnames(choicematrix) = altnames
      cat('Overview of choices for MNL model component:\n')
      print(round(choicematrix,2))
      cat("\n")
      if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
    }

    return(invisible(TRUE))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    # Store useful values
    nObs  <- length(choice)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choice))
    choice[!rows]=alternatives[1]

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

    # Store useful values
    nObs  <- length(choice)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choice))
    choice[!rows]=alternatives[1]

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
                             FUN=function(i) (choice==altcodes[i])*V[[altnames[i]]])
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
      P <- lapply(P, function(p) {
        # change likelihood of excluded columns for 1
        if(is.vector(p)) p[!rows]  <- NA
        if(is.matrix(p)) p[!rows,] <- NA
        if(is.array(p) & length(dim(p))==3) p[!rows,,] <- NA
        return(p)
      })
      # add an additional column with chosen
      P[["chosen"]]= 1/denom
      if(is.vector(P[["chosen"]])) P[["chosen"]][!rows]  <- NA
      if(is.matrix(P[["chosen"]])) P[["chosen"]][!rows,] <- NA
      if(is.array(P[["chosen"]]) && length(dim(P[["chosen"]]))==3) P[["chosen"]][!rows,,] <- NA
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
    nObs  <- length(choice)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choice))
    choice[!rows]=alternatives[1]

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
                             FUN=function(i) (choice==altcodes[i])*V[[altnames[i]]])
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
    availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
    choicematrix = matrix(0,nrow=4,ncol=length(altnames))
    choicematrix[1,] = availprint
    j=1
    while(j<= length(altnames)){
      choicematrix[2,j]=sum(choice==altcodes[j] & rows) # number of times each alt is chosen
      j=j+1
    }
    choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
    choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
    rownames(choicematrix) = c("Times available","Times chosen","Percentage of choice overall","Percentage of choice when available")
    colnames(choicematrix) = altnames

    # write "choicematrix" to a file named "modelName_tempOutput.txt"
    apollo_control <- get("apollo_control",parent.frame(), inherits=TRUE)
    fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    fileName <- file.path(tempdir(),fileName)
    tryCatch( fileConn <- file(fileName, open="at"), error=function(e) fileConn <- NA )
    if(!anyNA(fileConn)){
      sink(fileConn)
      if(apollo_control$noDiagnostics==FALSE){
        cat('Overview of choices for MNL model component:\n')
        print(round(choicematrix,2))
        cat("\n")}
      if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
      if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      sink()
      close(fileConn)
    } else cat('Could not write choice stats of MNL model in temporary file.\n')

    return(P)
  }
}
