#' Calculates probabilities of a cross nested logit
#'
#' Calculates probabilities of a cross nested logit model.
#'
#' For the model to be consistent with utility maximisation, the estimated value of the lambda parameter of all nests
#' should be between 0 and 1. Lambda parameters are inversely proportional to the correlation between the error terms of 
#' alternatives in a nest. If lambda=1,  there is no relevant correlation between the unobserved
#' utility of alternatives in that nest.
#' The tree must contain an upper nest called \code{"root"}. The lambda parameter of the root is automatically
#' set to 1 if not specified in \code{nlNests}. And while setting it to another value is possible, it is not
#' recommended.
#' Alpha parameters inside \code{cnlStructure} should be between 0 and 1. Using a transformation to ensure
#' this constraint is satisfied is recommended (e.g. logistic transformation).
#' @param cnl_settings List of inputs of the CNL model. It should contain the following.
#'                     \itemize{
#'                       \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item choiceVar: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item V: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item cnlNests: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named according to nests. The lambda at the root is fixed to 1, and therefore does not need to be defined.
#'                       \item cnlStructure: Numeric matrix. One row per nest and one column per alternative. Each element of the matrix is the alpha parameter of that (nest, alternative) pair.
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
#' @importFrom stats setNames
#' @examples
#' ### Load data
#' data(apollo_modeChoiceData)
#' database <- apollo_modeChoiceData
#' rm(apollo_modeChoiceData)
#'
#' ### Parameters
#' b = list(asc_1=0, asc_2=0, asc_3=0, asc_4=0, tt=0, tc=0, acc=0,
#'          lambda_fastPT=0.5, lambda_groundPT=0.5, alpha_rail_fastPT=0.5)
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
#' cnlStructure     = matrix(0, nrow=3, ncol=4)
#' cnlStructure[1,] = c( 0,  0, 1, b$alpha_rail_fastPT  ) # fastPT
#' cnlStructure[2,] = c( 0,  1, 0, 1-b$alpha_rail_fastPT) # groundPT
#' cnlStructure[3,] = c( 1,  0, 0, 0                    ) # car
#'
#' ### CNL settings
#' cnlSettings <- list(
#'    alternatives = c(car=1, bus=2, air=3, rail=4),
#'    avail        = list(car=database$av_car, bus=database$av_bus,
#'                        air=database$av_air, rail=database$av_rail),
#'    choiceVar    = database$choice,
#'    V            = V,
#'    cnlNests     = list(fastPT=b$lambda_fastPT, groundPT=b$lambda_groundPT, car=1),
#'    cnlStructure = cnlStructure
#' )
#'
#' ### Compute choice probabilities using CNL model
#' apollo_cnl(cnlSettings, functionality="estimate")
#' @export
apollo_cnl <- function(cnl_settings, functionality){
  if(is.null(cnl_settings[["alternatives"]])) stop("The cnl_settings list needs to include an object called \"alternatives\"!")
  if(is.null(cnl_settings[["avail"]])) stop("The cnl_settings list needs to include an object called \"avail\"!")
  if(is.null(cnl_settings[["choiceVar"]])) stop("The cnl_settings list needs to include an object called \"choiceVar\"!")
  if(is.null(cnl_settings[["V"]])) stop("The cnl_settings list needs to include an object called \"V\"!")
  if(is.null(cnl_settings[["cnlNests"]])) stop("The cnl_settings list needs to include an object called \"cnlNests\"!")
  if(is.null(cnl_settings[["cnlStructure"]])) stop("The cnl_settings list needs to include an object called \"cnlStructure\"!")
  if(is.null(cnl_settings[["rows"]])) cnl_settings[["rows"]]="all"

  alternatives = cnl_settings[["alternatives"]]
  avail        = cnl_settings[["avail"]]
  choiceVar    = cnl_settings[["choiceVar"]]
  V            = cnl_settings[["V"]]
  cnlNests     = cnl_settings[["cnlNests"]]
  cnlStructure = cnl_settings[["cnlStructure"]]
  rows         = cnl_settings[["rows"]]

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    # Store useful values
    apollo_control <- tryCatch(get("apollo_control", parent.frame(), inherits=FALSE),
                               error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))

    nObs  <- length(choiceVar)
    nAlts <- length(V)
    nameAltV <- names(V)
    avail_set <- FALSE
    altnames <- names(alternatives)
    altcodes <- alternatives
    nestnames<- names(cnlNests)
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
      #warning("Full availability of alternatives assumed for the CNL component.")
    }

    # Reorder V to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]

    if(apollo_control$noValidation==FALSE){
      # Checks that user doesn't attempt to change the root lambda
      if("root" %in% names(cnlNests)) stop("The root should not be included in argument cnlNests.")

      # Check there are at least three alternatives
      if(nAlts<3) stop("CNL requires at least three alternatives")

      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")

      # Check that labels in choice match those in the utilities and availabilities
      choiceLabs <- unique(choiceVar)
      if(!all(altnames %in% names(V))) stop("Alternative labels in \"altnames\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternative labels in \"altnames\" do not match those in \"avail\".")

      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop("Value in choice column that is not included in altcodes.")

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

      # checks that are specific to cnlStructure component
      if(nrow(cnlStructure)!=length(nestnames)) stop("Tree structure needs one row per nest!")
      if(ncol(cnlStructure)!=nAlts) stop("Tree structure needs one column per alternative!")
      if(max(colSums(cnlStructure))>1) stop("Allocation parameters for some alternatives sum to value above 1!")
      if(min(colSums(cnlStructure))>1) stop("Allocation parameters for some alternatives sum to value below 1!")

      # confirm all checks are passed, print output and return TRUE
      #cat("\nAll checks passed for CNL model component\n")

    }

    if(apollo_control$noDiagnostics==FALSE){

      #if(avail_set==TRUE) warning("Availability not provided to 'apollo_cnl' (or some elements are NA). Full availability assumed.")

      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)

      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) 
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows) 
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = altnames
      #cat('Overview of choices for CNL model component:\n')
      #print(round(choicematrix,2))
      #cat("\n")
      #if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      #if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      #
      ## Print structure of the CNL model only if lambda is common for all sample.
      #if(all(sapply(cnlNests, function(x) length(x)==1))){
      #  out_tree = cbind(cnlStructure, unlist(cnlNests))
      #  out_tree = apply(out_tree, MARGIN=2, function(x){
      #    if(all(x %in% 0:1)) round(x,0) else round(x,4)
      #    return(x)
      #  } )
      #  rownames(out_tree)=nestnames
      #  colnames(out_tree)=c(altnames,"lambda")
      #  
      #  cat('Initial structure for CNL model component:\n')
      #  print(out_tree)
      #  cat("\n")
      #}
      
      content <- list(round(choicematrix,2),
                      ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
                      ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", "") )
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      if( all(sapply(cnlNests, function(x) length(x)==1)) ){
        out_tree = cbind(cnlStructure, unlist(cnlNests))
        out_tree = apply(out_tree, MARGIN=2, function(x){
          if(all(x %in% 0:1)) round(x,0) else round(x,4)
          return(x)
        } )
        rownames(out_tree)=nestnames
        colnames(out_tree)=c(altnames,"lambda")
        content[[length(content) + 1]] <- "Initial structure for CNL model component:"
        content[[length(content) + 1]] <- out_tree
      }
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for CNL model component:", content, apolloLog)
      
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
    altnames <- names(alternatives)
    altcodes <- alternatives
    nestnames<- names(cnlNests)
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
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]
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
    if(functionality=="raw" && length(choiceVar)==1 && is.na(choiceVar)) choiceVar = alternatives[1]
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames <- names(alternatives)
    altcodes <- alternatives
    nestnames<- names(cnlNests)
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
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]


    # extract chosen utility
    chosenV <- Reduce('+', lapply(as.list(1:nAlts),
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

    # work out denominator for within nest probs

    denom_within=list()
    t=1
    nests=nrow(cnlStructure)

    while(t<=nests){
      denom_within[[t]]=0
      j=1
      while(j<=nAlts){
        denom_within[[t]]=denom_within[[t]]+(cnlStructure[t,j]*eV[[altnames[j]]])^(1/cnlNests[[t]])
        j=j+1
      }
      t=t+1
    }

    # work out within-nest probs

    Pwithin=list()

    j=1
    while(j<=nAlts){
      Pwithin[[j]]=list()
      t=1
      while(t<=nests){
        Pwithin[[j]][[t]]=(cnlStructure[t,j]*eV[[altnames[j]]])^(1/cnlNests[[t]])/(denom_within[[t]]+(denom_within[[t]]==0)) # includes a failsafe in denominator for empty nests, numerator will ensure ratio is 0 anyway
        t=t+1
      }
      j=j+1
    }

    # work out nest probs

    Pnest=list()
    denom_nest=0
    t=1
    while(t<=nests){
      denom_nest=denom_nest+denom_within[[t]]^cnlNests[[t]]
      t=t+1}
    t=1
    while(t<=nests){
      Pnest[[t]]=(denom_within[[t]]^cnlNests[[t]])/denom_nest
      t=t+1}

    # work individual probs

    Palts=list()
    j=1
    while(j<=nAlts){
      Palts[[j]]=0
      t=1
      while(t<=nests){
        Palts[[j]]=Palts[[j]]+Pnest[[t]]*Pwithin[[j]][[t]]
        t=t+1}
      j=j+1
    }

    names(Palts)=names(V)

    if(functionality=="prediction"|(functionality=="raw")){
      P=Palts
      # add an additional column with chosen
      chosenP <- (choiceVar==altcodes[1])*Palts[[altnames[1]]]
      for(i in 2:nAlts) chosenP <- chosenP + (choiceVar==altcodes[i])*Palts[[altnames[i]]]
      if(functionality=="prediction") P[["chosen"]]=chosenP
      P <- lapply(P, function(p) {
        # change likelihood of excluded columns for 1
        if(is.vector(p)) p[!rows]  <- ifelse(functionality=="prediction",NA,1)
        if(is.matrix(p)) p[!rows,] <- ifelse(functionality=="prediction",NA,1)
        if(is.array(p) & length(dim(p))==3) p[!rows,,] <- ifelse(functionality=="prediction",NA,1)
        return(p)
      })
    } else {
      # keep only probability for chosen alternative
      chosenP <- (choiceVar==altcodes[1])*Palts[[altnames[1]]]
      for(i in 2:nAlts) chosenP <- chosenP + (choiceVar==altcodes[i])*Palts[[altnames[i]]]
      P=chosenP
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

    P = apollo_cnl(cnl_settings, functionality="estimate")
    
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames <- names(alternatives)
    altcodes <- alternatives
    nestnames<- names(cnlNests)
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), names(alternatives))
    }
    
    ### Reorder V and avail to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]


    if(avail_set) availprint = rep(length(choiceVar),length(altnames)) else{
      availprint=colSums(rows*matrix(unlist(avail), ncol = length(avail)))
    }

    choicematrix=matrix(0,nrow=4,ncol=length(altnames))

    choicematrix[1,]=availprint
    j=1
    while(j<= length(altnames)){
      choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows)
      j=j+1
    }

    choicematrix[3,] = choicematrix[2,]/sum(rows)*100
    choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100
    choicematrix[4,!is.finite(choicematrix[4,])] <- 0
    rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
    colnames(choicematrix) = altnames
    
    # Print structure of the CNL model only if lambda is common for all sample.
    #if(all(sapply(cnlNests, function(x) length(x)==1))){
    #  out_tree = cbind(cnlStructure, unlist(cnlNests))
    #  out_tree = apply(out_tree, MARGIN=2, function(x){
    #    if(all(x %in% 0:1)) round(x,0) else round(x,4)
    #    return(x)
    #  } )
    #  rownames(out_tree)=nestnames
    #  colnames(out_tree)=c(altnames,"lambda")
    #  
    #  # write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
    #  apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
    #                              error=function(e){
    #                                cat("apollo_cnl could not retrieve apollo_control. No diagnostics in output.\n")
    #                                return(NA)
    #                              } )
    #  if(!(length(apollo_control)==1 && is.na(apollo_control))){
    #    fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    #    fileName <- file.path(tempdir(),fileName)
    #    fileConn <- tryCatch( file(fileName, open="at"),
    #                          error=function(e){
    #                            cat('apollo_cnl could not write diagnostics to temporary file. No diagnostics in output.\n')
    #                            return(NA)
    #                          })
    #    if(!anyNA(fileConn)){
    #      sink(fileConn)
    #      on.exit({if(sink.number()>0) sink(); close(fileConn)})
    #      if(apollo_control$noDiagnostics==FALSE){
    #        cat('Overview of choices for CNL model component:\n')
    #        print(round(choicematrix,0))
    #        cat("\n")
    #      }
    #      cat("Final structure for CNL model component:\n")
    #      print(out_tree)
    #      cat("\n")
    #      if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
    #    }
    #  }
    #}
    
    content <- list(round(choicematrix,2))
    if(any(choicematrix[4,]==0)) content[[length(content)+1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(choicematrix[4,]==1)) content[[length(content)+1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    if( all(sapply(cnlNests, function(x) length(x)==1)) ){
      out_tree = cbind(cnlStructure, unlist(cnlNests))
      out_tree = apply(out_tree, MARGIN=2, function(x){
        if(all(x %in% 0:1)) round(x,0) else round(x,4)
        return(x)
      } )
      rownames(out_tree)=nestnames
      colnames(out_tree)=c(altnames,"lambda")
      content[[length(content) + 1]] <- "Final structure for CNL model component:"
      content[[length(content) + 1]] <- out_tree
    }
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for CNL model component:", content, apolloLog)
    
    
    return(P)
  }
}
