#' Calculates Multinomial Logit probabilities
#'
#' Calculates the probabilities of a Multinomial Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' @param mnl_settings List of inputs of the MNL model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}. Set to \code{"all"} by default if omitted.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                     }
#' @param functionality Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call \code{apollo_probabilities}, though the user can also call \code{apollo_probabilities} manually with a given functionality for testing/debugging. Possible values are:
#'                      \itemize{
#'                        \item \strong{\code{"components"}}: For further processing/debugging, produces likelihood for each model component (if multiple components are present), at the level of individual draws and observations.
#'                        \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
#'                        \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"gradient"}}: For model estimation, produces analytical gradients of the likelihood, where possible.
#'                        \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
#'                        \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
#'                        \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
#'                        \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
#'                        \item \strong{\code{"report"}}: Prepares output summarising model and choiceset structure.
#'                        \item \strong{\code{"shares_LL"}}: Produces overall model likelihood with constants only.
#'                        \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"zero_LL"}}: Produces overall model likelihood with all parameters at zero.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"gradient"}}: List containing the likelihood and gradient of the model component.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{mnl_settings}.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Choice overview
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @export
#' @importFrom utils capture.output
apollo_mnl <- function(mnl_settings, functionality){
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE, silent=FALSE, analyticGrad=TRUE)) ))
  
  ### Set or extract componentName
  modelType   = "MNL"
  if(is.null(mnl_settings[["componentName"]])){
    mnl_settings[["componentName"]] = ifelse(!is.null(mnl_settings[['componentName2']]),
                                             mnl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && mnl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', mnl_settings[["componentName"]], '" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, mnl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", mnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used
  if(!is.null(mnl_settings[["utilities"]])) names(mnl_settings)[which(names(mnl_settings)=="utilities")]="V"
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  if( !is.null(apollo_inputs[[paste0(mnl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load mnl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(mnl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded mnl_settings, restore the one received as argument
    if(is.null(tmp$V)) tmp$V <- mnl_settings$V
    mnl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    if(is.null(names(mnl_settings))) stop('SYNTAX ISSUE - All elements inside the inputs lists for model components must be named, e.g. mnl_settings=list(alternatives=c(...), avail=...).')
    if(is.null(mnl_settings[["componentName"]])) stop('SYNTAX ISSUE - The settings of at least one model component is missing the mandatory "componentName" object.')
    
    # functionality
    test <- functionality %in% c("estimate","prediction","validate","zero_LL","shares_LL","conditionals","output","raw","preprocess", "components", "gradient","hessian", "report")
    if(!test) stop("SYNTAX ISSUE - Non-permissable setting for \"functionality\" for model component \"",mnl_settings$componentName,"\"")
    
    # Check for mandatory inputs
    mandatory <- c("alternatives", "choiceVar", "V")
    for(i in mandatory) if(!(i %in% names(mnl_settings))) stop('SYNTAX ISSUE - The inputs list for model component "', mnl_settings$componentName, 
                                                               '" needs to include an object called "', i,'"!')
    # Check for optional inputs (avail and rows)
    if(is.null(mnl_settings[["rows"]])) mnl_settings[["rows"]]="all"
    if(is.null(mnl_settings[['avail']])){
      mnl_settings[['avail']]=1
      if(!apollo_inputs$silent && functionality=='validate') apollo_print('Setting "avail" is missing, so full availability is assumed.')
    }
    
    ### Store useful values
    mnl_settings$altnames = names(mnl_settings$alternatives)
    mnl_settings$altcodes = mnl_settings$alternatives
    mnl_settings$nAlt     = length(mnl_settings$alternatives)
    mnl_settings$nObs <- tryCatch(if(!is.null(apollo_inputs$database)) nrow(apollo_inputs$database) else stop('x'),
                                  error=function(e){
                                    lenV <- sapply(mnl_settings$V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                                    lenA <- sapply(mnl_settings$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                                    lenC <- ifelse(!is.null(mnl_settings$choiceVar),length(mnl_settings$choiceVar),length(mnl_settings$choiceShares[[1]]))
                                    return(max(lenV, lenA, lenC))
                                  })
    
    ### Format checks
    # alternatives
    test <- is.vector(mnl_settings$alternatives) & !is.null(names(mnl_settings$alternatives))
    if(!test) stop("SYNTAX ISSUE - The \"alternatives\" argument for model component \"",mnl_settings$componentName,"\" needs to be a named vector")
    # avail
    test <- is.list(mnl_settings$avail) || (length(mnl_settings$avail)==1 && mnl_settings$avail==1)
    if(!test) stop("SYNTAX ISSUE - The \"avail\" argument for model component \"",mnl_settings$componentName,"\" needs to be a list or set to 1")
    if(is.list(mnl_settings$avail)){
      lenA <- sapply(mnl_settings$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
      test <- all(lenA==mnl_settings$nObs | lenA==1)
      if(!test) stop("SYNTAX ISSUE - All entries in \"avail\" for model component \"",mnl_settings$componentName,"\" need to be a scalar or a vector with one entry per observation in the \"database\"")
    }
    test <- is.vector(mnl_settings$choiceVar) && (length(mnl_settings$choiceVar)==mnl_settings$nObs || length(mnl_settings$choiceVar)==1)
    if(!test) stop("SYNTAX ISSUE - The \"choiceVar\" argument for model component \"",mnl_settings$componentName,"\" needs to be a scalar or a vector with one entry per observation in the \"database\"")
    
    # rows
    test <- is.vector(mnl_settings$rows)
    #test <- test && ( (is.logical(mnl_settings$rows) && length(mnl_settings$rows)==mnl_settings$nObs) || (length(mnl_settings$rows)==1 && mnl_settings$rows=="all") )
    #if(!test) stop("SYNTAX ISSUE - The \"rows\" argument for model component \"",mnl_settings$componentName,"\" needs to be \"all\" or a vector of logical statements with one entry per observation in the \"database\"")
    test <- test && ( (is.logical(mnl_settings$rows) || all(mnl_settings$rows%in%c(0,1))) && length(mnl_settings$rows)==mnl_settings$nObs) || ( (all(mnl_settings$rows%in%c(0,1)) && length(mnl_settings$rows)==mnl_settings$nObs) || (length(mnl_settings$rows)==1 && mnl_settings$rows=="all") )
    if(!test) stop("SYNTAX ISSUE - The \"rows\" argument for model component \"",mnl_settings$componentName,"\" needs to be \"all\" or a vector of logical statements or 0/1 entries with one entry per observation in the \"database\"")
    if( all(mnl_settings$rows %in% c(0,1)) ) (mnl_settings$rows <- mnl_settings$rows>0)
    
    ### Expand availabilities if necessary
    mnl_settings$avail_set <- FALSE
    if(length(mnl_settings$avail)==1 && mnl_settings$avail==1){
      mnl_settings$avail <- as.list(setNames(rep(1,mnl_settings$nAlt), mnl_settings$altnames))
      mnl_settings$avail_set <- TRUE
    }
    
    ### Check that avail and V are available for all alternatives
    if(!all(mnl_settings$altnames %in% names(mnl_settings$V))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"", mnl_settings$componentName,"\" do not match those in \"utilities\".")
    if(!all(mnl_settings$altnames %in% names(mnl_settings$avail))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"", mnl_settings$componentName,"\" do not match those in \"avail\".")
    if(length(mnl_settings$V)>length(mnl_settings$altnames)) stop("SYNTAX ISSUE - More utilities have been defined than there are alternatives for model component \"", mnl_settings$componentName,"\"!")
    if(length(mnl_settings$avail)>length(mnl_settings$altnames)) stop("SYNTAX ISSUE - More availabilities have been defined than there are alternatives for model component \"", mnl_settings$componentName,"\"!")
      
    ### Reorder availabilities and V
    mnl_settings$avail <- mnl_settings$avail[mnl_settings$altnames]
    mnl_settings$V     <- mnl_settings$V[mnl_settings$altnames]
    
    ### Expand rows if necessary, and update nObs
    if(length(mnl_settings$rows)==1 && mnl_settings$rows=="all") mnl_settings$rows <- rep(TRUE, mnl_settings$nObs)
    mnl_settings$nObs <- sum(mnl_settings$rows)
    # Filter rows, except for V
    if(any(!mnl_settings$rows)){
      mnl_settings$avail <- lapply(mnl_settings$avail, 
                                   function(av) if(length(av)==1) return(av) else return(av[mnl_settings$rows]))
      mnl_settings$choiceVar <- apollo_keepRows(mnl_settings$choiceVar, mnl_settings$rows)
    }
    
    ### Create Y
    mnl_settings$Y <- lapply(as.list(mnl_settings$alternatives), function(i) mnl_settings$choiceVar==i)
    
    # Record availability of chosen alternative
    mnl_settings$chosenAvail <- Reduce('+', mapply('*', mnl_settings$Y, mnl_settings$avail, SIMPLIFY=FALSE))
    
    # Determine which mnl likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation for MNL available")
    # Using R likelihood
    mnl_settings$probs_MNL=function(mnl_settings, all=FALSE){
      # Fix choiceVar if "raw" and choiceVar==NA
      mnl_settings$choiceNA = FALSE
      if(all(is.na(mnl_settings$choiceVar))){
        mnl_settings$choiceVar = mnl_settings$alternatives[1]
        mnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      V <- mapply(function(v,a) apollo_setRows(v, !a, 0), mnl_settings$V, mnl_settings$avail, SIMPLIFY=FALSE)
      # if probabilities for all alternatives are requested, then P is a list
      if(all){
        if(apollo_inputs$apollo_control$subMaxV){
          ### work with subtracting the maxV
          maxV <- do.call(pmax, V)
          V <- lapply(V, "-", maxV)
          V <- lapply(X=V, FUN=exp)
          V <- mapply('*', V, mnl_settings$avail, SIMPLIFY = FALSE)
          denom = Reduce('+', V)
          P <- lapply(V, "/", denom)
        } else {
          ### work with subtracting the chosenV
          chosenV <- mapply("*", mnl_settings$Y, V, SIMPLIFY=FALSE)
          chosenV <- Reduce('+', chosenV)
          V <- lapply(X=V, "-", chosenV)
          V <- lapply(X=V, FUN=exp)
          # consider availabilities (it assumes V and avail are in the same order)
          V <- mapply('*', V, mnl_settings$avail, SIMPLIFY = FALSE)
          # calculate the denominator of the Logit probability expression
          denom = Reduce('+',V)
          P <- lapply(V, "/", denom)
          if(any(sapply(P, anyNA))){
            P <- lapply(P, function(p){p[is.na(p)] <- 1; return(p)} )
            P <- mapply('*', P, mnl_settings$avail, SIMPLIFY = FALSE)
            sP <- Reduce("+", P)
            P <- lapply(P, "/", sP)
          }
        }
        if(!mnl_settings$choiceNA) P[["chosen"]] <- Reduce("+", mapply("*", mnl_settings$Y, P, SIMPLIFY=FALSE))
      } else { 
        # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
        ### work with subtracting the chosenV
        chosenV <- mapply("*", mnl_settings$Y, V, SIMPLIFY=FALSE)
        chosenV <- Reduce('+', chosenV)
        V <- lapply(X=V, "-", chosenV)
        V <- lapply(X=V, FUN=exp)
        # consider availabilities (it assumes V and avail are in the same order)
        V <- mapply('*', V, mnl_settings$avail, SIMPLIFY = FALSE)
        # calculate the denominator of the Logit probability expression
        denom = Reduce('+', V)
        P <- mnl_settings$chosenAvail/denom
      }
      return(P)
    }
    
    # Create diagnostics function for MNL
    mnl_settings$mnl_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
      # turn scalar availabilities into vectors
      for(i in 1:length(inputs$avail)) if(length(inputs$avail[[i]])==1) inputs$avail[[i]] <- rep(inputs$avail[[i]], inputs$nObs)
      
      # Construct summary table of availabilities and market share
      choicematrix = matrix(0, nrow=4, ncol=length(inputs$altnames), 
                            dimnames=list(c("Times available", "Times chosen", "Percentage chosen overall",
                                            "Percentage chosen when available"), inputs$altnames))
      choicematrix[1,] = unlist(lapply(inputs$avail, sum))
      for(j in 1:length(inputs$altnames)) choicematrix[2,j] = sum(inputs$choiceVar==inputs$altcodes[j]) # number of times each alt is chosen
      choicematrix[3,] = choicematrix[2,]/inputs$nObs*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      
      if(!apollo_inputs$silent & data){
        if(any(choicematrix[4,]==0)) apollo_print("Some alternatives are never chosen in your data!", type="w")
        if(any(choicematrix[4,]>=100)) apollo_print("Some alternatives are always chosen when available!", type="w")
        #if(inputs$avail_set) apollo_print("Availability not provided (or some elements are NA). Full availability assumed.", type="i")
        apollo_print("\n")
        apollo_print(paste0('Overview of choices for MNL model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        print(round(choicematrix,2))
        cat("\n")
      }
      
      return(invisible(TRUE))
    }
    
    # Store model type
    mnl_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient", "hessian"))
    test <- test && all(sapply(mnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    mnl_settings$gradient <- FALSE
    if(test){
      mnl_settings$V        <- mnl_settings$V[mnl_settings$altnames] # reorder V
      mnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, mnl_settings$V)
      mnl_settings$gradient <- !is.null(mnl_settings$dV)
    }; rm(test)
    
    # Construct necessary input for hessian
    test <- !is.null(mnl_settings$gradient) && mnl_settings$gradient && apollo_inputs$apollo_control$analyticHessian
    mnl_settings$hessian <- test
    if(test){
      mnl_settings$ddV <- list()
      alts=mnl_settings$altnames
      pars=names(mnl_settings$dV)
      for(k1 in pars){
        mnl_settings$ddV[[k1]] <- list() 
        for(k2 in pars){
          mnl_settings$ddV[[k1]][[k2]] <- list() 
          for(j in alts){
            mnl_settings$ddV[[k1]][[k2]][[j]] <- Deriv::Deriv(f=mnl_settings$dV[[k1]][[j]], x=k2)
            if(is.null(mnl_settings$ddV[[k1]][[k2]][[j]])) mnl_settings$hessian <- FALSE
          }
        }
      }
    }; rm(test)
    
    # Return mnl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      mnl_settings$V      <- NULL
      return(mnl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(mnl_settings$V, is.function))){
    mnl_settings$V = lapply(mnl_settings$V, function(f) if(is.function(f)) f() else f )
  } 
  mnl_settings$V <- lapply(mnl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V and drop rows if necessary
  mnl_settings$V <- mnl_settings$V[mnl_settings$altnames]
  if(!all(mnl_settings$rows)) mnl_settings$V <- lapply(mnl_settings$V, apollo_keepRows, r=mnl_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Check that alternatives are named in altcodes and V
    if(is.null(mnl_settings$altnames) || is.null(mnl_settings$altcodes) || is.null(names(mnl_settings$V))) stop("SYNTAX ISSUE - Alternatives for model component \"",mnl_settings$componentName,"\" must be named, both in 'alternatives' and 'utilities'.")
    
    if(!apollo_inputs$apollo_control$noValidation){
      # Check there are no repeated alternatives names
      if(length(unique(mnl_settings$altnames))!=length(mnl_settings$altnames)) stop('SYNTAX ISSUE - Names of alternatives must be unique. Check definition of "alternatives".')
      
      # Check that there are at least two alternatives
      minAlts <- 2
      if(mnl_settings$nAlt<minAlts) stop("SYNTAX ISSUE - Model component \"",mnl_settings$componentName,"\"  requires at least ", minAlts, " alternatives")
      
      # Check if LLC is required with large numbers of alternatives
      if(mnl_settings$nAlt>20 && !is.null(apollo_inputs$apollo_control$calculateLLC) && apollo_inputs$apollo_control$calculateLLC) apollo_print(paste0("The number of the alternatives for model component \"",mnl_settings$componentName,"\" is large and you may consider setting apollo_control$calculateLLC=FALSE to avoid the calculation of the log-likelihod with constants only."), pause=0, type="i")

      # Check that choice vector is not empty
      if(length(mnl_settings$choiceVar)==0) stop("SYNTAX ISSUE - Choice vector is empty for model component \"",mnl_settings$componentName,"\"")
      
      if(mnl_settings$nObs==0) stop("INPUT ISSUE - No data for model component \"",mnl_settings$componentName,"\"")
      
      # Check V and avail elements are named correctly
      if(!all(mnl_settings$altnames %in% names(mnl_settings$V))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"",mnl_settings$componentName,"\" do not match those in \"utilities\".")
      if(!all(mnl_settings$altnames %in% names(mnl_settings$avail))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"",mnl_settings$componentName,"\" do not match those in \"avail\".")
      if(length(mnl_settings$V)>length(mnl_settings$altnames)) stop("SYNTAX ISSUE - More utilities have been defined than there are alternatives for model component \"", mnl_settings$componentName,"\"!")
      if(length(mnl_settings$avail)>length(mnl_settings$altnames)) stop("SYNTAX ISSUE - More availabilities have been defined than there are alternatives for model component \"", mnl_settings$componentName,"\"!")
      
      # Check that there are no values in the choice column for undefined alternatives
      mnl_settings$choiceLabs <- unique(mnl_settings$choiceVar)
      if(!all(mnl_settings$choiceLabs %in% mnl_settings$altcodes)) stop("INPUT ISSUE - The data contains values for \"choiceVar\" for model component \"",mnl_settings$componentName,"\" that are not included in \"alternatives\".")
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(mnl_settings$avail)) if( !all(unique(mnl_settings$avail[[i]]) %in% 0:1) ) stop("INPUT ISSUE - Some availability values for model component \"",mnl_settings$componentName,"\" are not 0 or 1.")
      # check that at least 2 alternatives are available in at least one observation
      if(max(Reduce('+',mnl_settings$avail))==1) stop("SPECIFICATION ISSUE - Only one alternative is available for each observation for model component \"",mnl_settings$componentName,"!")
      # check that nothing unavailable is chosen
      for(j in 1:mnl_settings$nAlt) if(any(mnl_settings$choiceVar==mnl_settings$altcodes[j] & mnl_settings$avail[[j]]==0)){
        #stop("SPECIFICATION ISSUE - The data contains cases where alternative ",
        #     mnl_settings$altnames[j]," is chosen for model component \"",
        #     mnl_settings$componentName,"\" despite being listed as unavailable\n")
        txt <-  paste0('The data contains cases where alternative ', mnl_settings$altnames[j], 
                       ' is chosen for model component "',mnl_settings$componentName, '" despite being', 
                       ' listed as unavailable. This will cause the chosen probability to be', 
                       ' zero, and potentially lead to an invalid LL.')
        apollo_print(txt, type="w")
      }
      
      # Check that no available alternative has utility = NA
      # Requires setting non available alternatives utility to 0 first
      mnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), mnl_settings$V, mnl_settings$avail, SIMPLIFY=FALSE)
      if(!all(sapply(mnl_settings$V, function(v) all(is.finite(v))))) stop('CALCULATION ISSUE - Some utilities for model component "',
                                                                     mnl_settings$componentName, 
                                                                     '" contain values that are not finite numbers!')
    } 
    
    if(!apollo_inputs$apollo_control$noDiagnostics) mnl_settings$mnl_diagnostics(mnl_settings, apollo_inputs)
    
    testL = mnl_settings$probs_MNL(mnl_settings, all=FALSE)
    if(any(!mnl_settings$rows)) testL <- apollo_insertRows(testL, mnl_settings$rows, 1)
    if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",mnl_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",mnl_settings$componentName,"\"", sep=""))
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    for(i in 1:length(mnl_settings$avail)) if(length(mnl_settings$avail[[i]])==1) mnl_settings$avail[[i]] <- rep(mnl_settings$avail[[i]], mnl_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, mnl_settings$avail)) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    if(any(!mnl_settings$rows)) P <- apollo_insertRows(P, mnl_settings$rows, 1)
    return(P)
  }

  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality=="shares_LL"){
    for(i in 1:length(mnl_settings$avail)) if(length(mnl_settings$avail[[i]])==1) mnl_settings$avail[[i]] <- rep(mnl_settings$avail[[i]], mnl_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, mnl_settings$avail)) # number of available alts in each observation
    Y = do.call(cbind,mnl_settings$Y)
    if(var(nAvAlt)==0){
      Yshares = colSums(Y)/nrow(Y)
      P = as.vector(Y%*%Yshares)
    } else {
      ## Estimate model with constants only
      mnl_ll = function(b, A, Y) as.vector(Y%*%c(b,0) - log(rowSums( A%*%exp(c(b,0)) )))
      A = do.call(cbind, mnl_settings$avail)
      b = maxLik::maxLik(mnl_ll, start=rep(0, mnl_settings$nAlt - 1), 
                         method='BFGS', finalHessian=FALSE, A=A, Y=Y)$estimate
      P = exp(mnl_ll(b, A, Y))
    }
    if(any(!mnl_settings$rows)) P <- apollo_insertRows(P, mnl_settings$rows, 1)
    return(P)
  }
  
  # ################################################################# #
  #### functionality="estimate/prediction/conditionals/raw/output" ####
  # ################################################################# #
  
  if(functionality %in% c("estimate","conditionals", "components", "output")){
    P <- mnl_settings$probs_MNL(mnl_settings, all=FALSE)
    if(any(!mnl_settings$rows)) P <- apollo_insertRows(P, mnl_settings$rows, 1)
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- mnl_settings$probs_MNL(mnl_settings, all=TRUE)
    if(any(!mnl_settings$rows)) P <- lapply(P, apollo_insertRows, r=mnl_settings$r, val=NA)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(mnl_settings$dV)) stop("INTERNAL ISSUE - Analytical gradients cannot be calculated because the derivatives of the utilities are not available. Please set apollo_control$analyticGrad=FALSE.")
    for(k in 1:length(mnl_settings$dV)) if(!all( sapply(mnl_settings$dV[[k]], is.function) )) stop("INTERNAL ISSUE - Analytical gradients cannot be calculated because not al the derivatives of the utilities are functions. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_mnl could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_mnl could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities
    P    <- mnl_settings$probs_MNL(mnl_settings, all=TRUE)
    Pcho <- P[["chosen"]]
    P    <- P[-which(names(P)=="chosen")]
    
    # Calculate gradient
    J <- length(mnl_settings$dV[[1]]) # number of alternatives
    K <- length(mnl_settings$dV) # number of parameters
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    G <- setNames(vector(mode="list", length=K), names(mnl_settings$dV))
    r <- all(mnl_settings$rows) # TRUE if all rows are used (no rows excluded)
    a <- sapply(mnl_settings$avail, function(a) if(length(a)==1) a==1 else all(a==1)) # TRUE if all available
    for(k in 1:K){
      G[[k]] <- 0
      for(j in 1:J){
        # Calculate dVj/dbk, remove rows, expand it and replace unavailables rows by zero if necessary
        dVjk <- mnl_settings$dV[[k]][[j]]
        environment(dVjk) <- e
        dVjk <- dVjk()
        if(!r) dVjk <- apollo_keepRows(dVjk, mnl_settings$rows)
        if(length(dVjk)==1 && !a[j]) dVjk <- rep(dVjk, mnl_settings$nObs)
        if(!a[j]) dVjk <- apollo_setRows(dVjk, !mnl_settings$avail[[j]], 0)
        # Calculate gradient
        G[[k]] <- G[[k]] + (mnl_settings$Y[[j]] - P[[j]])*dVjk
      }
      G[[k]] <- Pcho*G[[k]]
      if(is.array(G[[k]])) rownames(G[[k]]) <- NULL else names(G[[k]]) <- NULL
    }; rm(dVjk)
    
    # Restore rows 
    if(!all(mnl_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, mnl_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=mnl_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############################# #
  #### functionality="hessian" ####
  # ############################# #
  if(functionality=="hessian"){
    ### TO DO
      # Checks before running
    
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - ",
                                                   "apollo_mnl could not ",
                                                   "fetch apollo_beta for ",
                                                   "hessian estimation."))
    
    ### Calculate necessary input
    J <- length(mnl_settings$dV[[1]]) # number of alternatives
    K <- length(mnl_settings$dV) # number of parameters
    N <- mnl_settings$nObs  # number of obs
    P <- mnl_settings$probs_MNL(mnl_settings, all=TRUE)
    L <- P[["chosen"]]
    L[L==0] <- 1e-50 # Remove zeros
    P <- P[-which(names(P)=="chosen")]
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, 
                    list(apollo_inputs=apollo_inputs)), hash=TRUE)
    r <- all(mnl_settings$rows) # TRUE if all rows are used (no rows excluded)
    A <- mnl_settings$avail
    a <- sapply(A, function(a) if(length(a)==1) a==1 else all(a==1)) # TRUE if all available
    pars <- names(mnl_settings$dV)
    alts <- names(mnl_settings$dV[[1]])
    
    ### Calculate gradient of probabilities for all alternatives and for chosen
    d1P <- list() ## derivatives of all P
    d1L <- setNames(vector(mode="list", length=K), pars) ## derivatives of chosen P
    for(k in 1:K){
      d1P[[k]] <- setNames(vector(mode="list", length=J), alts)
      for(j in 1:J) d1P[[k]][[j]]=0
      d1L[[k]] <- 0
      for(j in 1:J){
        # Calculate dVj/dbk, remove rows, expand it and replace unavailables rows by zero if necessary
        dVjk <- mnl_settings$dV[[k]][[j]]
        environment(dVjk) <- e
        dVjk <- dVjk()
        if(!r) dVjk <- apollo_keepRows(dVjk, mnl_settings$rows)
        if(length(dVjk)==1 && !a[j]) dVjk <- rep(dVjk, mnl_settings$nObs)
        if(!a[j]) dVjk <- apollo_setRows(dVjk, !mnl_settings$avail[[j]], 0)
        # calculate gradient of P
        for(i in 1:J){
          d1P[[k]][[i]] = d1P[[k]][[i]] + (ifelse(i==j,1,0) - P[[j]])*dVjk
        }
        # calculate gradient of L
        d1L[[k]] <- d1L[[k]] + (mnl_settings$Y[[j]] - P[[j]])*dVjk
      }
      d1L[[k]] <- L*d1L[[k]]
      d1P[[k]] <- mapply("*",d1P[[k]],P,SIMPLIFY=FALSE)
    }; rm(dVjk)
    
    # Calculate hessian of probability of chosen alternative
    d2L <- setNames(vector(mode="list", length=K), pars)
    for(k1 in 1:K){
      d2L[[k1]] <- setNames(vector(mode="list", length=K), pars)
      for(k2 in 1:k1){
        d2L[[k1]][[k2]] <- 0
        for(j in 1:J){
          yj  <- mnl_settings$Y[[j]]
          d1V <- mnl_settings$dV[[k1]][[j]]
          d2V <- mnl_settings$ddV[[k1]][[k2]][[j]]
          environment(d1V) <- e; environment(d2V) <- e
          d1V <- d1V()
          d2V <- d2V()
          if(!r){
            d1V <- apollo_keepRows(d1V, mnl_settings$rows)
            d2V <- apollo_keepRows(d2V, mnl_settings$rows)
          } 
          if(!a[j]){
            d1V <- apollo_setRows(d1V, !mnl_settings$avail[[j]], 0)
            d2V <- apollo_setRows(d2V, !mnl_settings$avail[[j]], 0)
          }
          # Update d2L only if d1V and d2V are not both zero.
          test <- is.vector(d1V) && length(d1V)==1 && d1V==0 &&
            is.vector(d2V) && length(d2V)==1 && d2V==0
          if(!test) d2L[[k1]][[k2]] <- d2L[[k1]][[k2]] + 
            (P[[j]] - yj)*d2V + d1V*d1P[[k2]][[j]]
        }
        d2L[[k1]][[k2]] <- d1L[[k2]]*d1L[[k1]]/L - L*d2L[[k1]][[k2]]
        # Restore rows
        if(!all(mnl_settings$rows)) d2L[[k1]][[k2]] <- apollo_insertRows(d2L[[k1]][[k2]], mnl_settings$rows, 0)
        d2L[[k2]][[k1]] <- d2L[[k1]][[k2]]
      }; rm(d1V, d2V)
    }; rm(d1P)
    
    # Restore rows in L and d1L (d2L already done above)
    if(!all(mnl_settings$rows)){
      L <- apollo_insertRows(L, mnl_settings$rows, 1)
      d1L <- lapply(d1L, apollo_insertRows, r=mnl_settings$rows, val=0)
    }
    
    # Return list with everything calculated
    return(list(like = L, grad=d1L, hess=d2L))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(mnl_settings$mnl_diagnostics(mnl_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(mnl_settings$mnl_diagnostics(mnl_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
