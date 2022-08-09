#' Calculates Exploded Logit probabilities
#'
#' Calculates the probabilities of an Exploded Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' 
#' The function calculates the probability of a ranking as a product of Multinomial Logit models with gradually reducing availability, where scale differences can be allowed for.
#' @param el_settings List of inputs of the Exploded Logit model. It shoud contain the following.
#'                    \itemize{
#'                     \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                     \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                     \item \strong{\code{choiceVars}}: List of numeric vectors. Contain choices for each position of the ranking. The list must be ordered with the best choice first, second best second, etc. It will usually be a list of columns from the database. Use value -1 if a stage does not apply for a given observations (e.g. when some individuals have shorter rankings).
#'                     \item \strong{\code{componentName}}: Character. Name given to model component.
#'                     \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                     \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                     \item \strong{\code{scales}}: List of vectors. Scale factors of each Logit model. At least one element should be normalized to 1. If omitted, scale=1 for all positions is assumed.
#'                    }
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
#'           \item \strong{\code{"prediction"}}: Not applicable (\code{NA}).
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{el_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"report"}}: Choice overview across stages.
#'           \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @importFrom stats setNames
#' @importFrom matrixStats rowCounts
#' @importFrom utils capture.output
#' @export
apollo_el <- function(el_settings, functionality){
  ### Set or extract componentName
  modelType   = "EL"
  if(is.null(el_settings[["componentName"]])){
    el_settings[["componentName"]] = ifelse(!is.null(el_settings[['componentName2']]),
                                            el_settings[['componentName2']], modelType)
    test <- functionality=="validate" && el_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 el_settings[['componentName']],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, el_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", el_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used
  if(!is.null(el_settings[["utilities"]])) names(el_settings)[which(names(el_settings)=="utilities")]="V"

  # ############################################### #
  #### Load pre-processing or do it if necessary ####
  # ############################################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(el_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    
    # Load el_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(el_settings$componentName, "_settings")]]
    # If there is no V inside the loaded el_settings, restore the one received as argument
    if(is.null(tmp$V)) tmp$V <- el_settings$V
    if(is.null(tmp$scales)) tmp$scales <- el_settings$scales
    el_settings <- tmp
    rm(tmp)
    
  } else { ### Do pre-processing
    ### Do pre-processing
    # Do pre-processing common to most models
    el_settings <- apollo_preprocess(inputs = el_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp && !apollo_inputs$silent) apollo_print("No C++ optimisation available for EL components.")
    # Using R likelihood
    el_settings$probs_EL <- function(el_settings){
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      el_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), el_settings$V, el_settings$avail[[1]], SIMPLIFY=FALSE)
      # Loop over stages, calculating the log-likelihood for each of them
      for(s in 1:el_settings$stages){
        # scale V's
        Vs <- lapply(el_settings$V, "*", el_settings$scales[[s]])
        # Substract V of chosen alternative to all other Vs and take their exponential
        Vi <- Reduce("+", mapply("*", el_settings$Y[[s]], Vs, SIMPLIFY=FALSE))
        Vs <- lapply(Vs, "-", Vi)
        Vs <- lapply(Vs, exp)
        # consider availabilities (it assumes V and avail are in the same order)
        Vs <- mapply('*', Vs, el_settings$avail[[s]], SIMPLIFY=FALSE)
        # calculate the denominator of the Logit probability expression
        denom <- Reduce('+', Vs)
        denom[el_settings$choiceVars[[s]]==-1 | !el_settings$chosenAvail[[s]]] <- 1
        if(s==1) P <- -log(denom) else P <- P - log(denom)
      }
      # Transform log-likelihood to likelihood
      P <- exp(P)
      return(P)
    }
    
    el_settings$el_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
      ### changes 28 July
        # Initialise summary table of availabilities and market share
        #choicematrix <- array(0, dim=c(4, inputs$nAlt+1, inputs$stages), 
        #                      dimnames=list(c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available"),
        #                                    c(inputs$altnames, "No choice"), paste("stage", 1:inputs$stages)) )
        choicematrix <- array(0, dim=c(4, inputs$nAlt, inputs$stages), 
                              dimnames=list(c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available"),
                                            c(inputs$altnames), paste("stage", 1:inputs$stages)) )
        # Calculate summary table for each stage and print it
        for(s in 1:inputs$stages){
          for(a in 1:inputs$nAlt){
            choicematrix[1,a,s] <- ifelse(length(inputs$avail[[s]][[a]])==1 && inputs$avail[[s]][[a]]==1, 
                                          inputs$nObs, sum(inputs$avail[[s]][[a]]) )
            choicematrix[2,a,s] <- sum(inputs$Y[[s]][[a]])
            choicematrix[3,a,s] <- choicematrix[2,a,s]/inputs$nObs*100
            choicematrix[4,a,s] <- choicematrix[2,a,s]/choicematrix[1,a,s]*100
            if(!is.finite(choicematrix[4,a,s])) choicematrix[4,a,s] <- 0
          }
          if(!apollo_inputs$silent & data){
            apollo_print("\n")
            apollo_print(paste0('Overview of choices for ', toupper(inputs$modeltype), ' model component', 
                                ifelse(inputs$componentName=='model', '', inputs$componentName), ', stage ', s,':'))
            print(round(choicematrix[,,s],2))
          }
        }
        
        if(!apollo_inputs$silent & data) for(a in 1:inputs$nAlt){
          if(sum(choicematrix[4,a,])==0) apollo_print(paste0('WARNING: Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'))
          if(choicematrix[4,a,1]==1) apollo_print(paste0('WARNING: Alternative "', inputs$altnames[a], '" is always chosen when available in model component "', inputs$componentName, '".'))
        }
        #if(inputs$avail_set==TRUE & !apollo_inputs$silent & data) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'))
      
        return(invisible(TRUE))
    }
    
    
    # Store model type
    el_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(el_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    el_settings$gradient <- FALSE
    if(test){
      el_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, el_settings$V)
      #el_settings$gradient <- !is.null(el_settings$dV)
    }; rm(test)
    
    # Return settings without V if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      el_settings$V      <- NULL
      if(!el_settings$fixedScales) el_settings$scales <- NULL
      return(el_settings)
    }
    
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  ### changes 28 July: this had mnl instead of el
  if(any(sapply(el_settings$V, is.function))) el_settings$V = lapply(el_settings$V, function(f) if(is.function(f)) f() else f)
  if(any(sapply(el_settings$scales, is.function))) el_settings$scales = lapply(el_settings$scales, function(f) if(is.function(f)) f() else f)
  el_settings$V <- lapply(el_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V and drop rows if neccesary
  el_settings$V <- el_settings$V[el_settings$altnames]
  if(!all(el_settings$rows)) el_settings$V <- lapply(el_settings$V, apollo_keepRows, r=el_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(el_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) el_settings$el_diagnostics(el_settings, apollo_inputs)
    
    testL <- el_settings$probs_EL(el_settings)
    if(any(!el_settings$rows)) testL <- apollo_insertRows(testL, el_settings$rows, 1)
    if(all(testL==0)) stop("All observations have zero probability at starting value for model component \"", el_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"", el_settings$componentName,"\"", sep=""))
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    P <- rep(1, el_settings$nObs)
    for(s in 1:el_settings$stages){
      nAvail <- Reduce("+", el_settings$avail[[s]])
      if(length(nAvail)>1) nAvail[el_settings$choiceVars[[s]]==-1 | nAvail==0] <- 1
      P <- P*1/nAvail
    }
    if(any(!el_settings$rows)) P <- apollo_insertRows(P, el_settings$rows, 1)
    return(P)
  }
  
  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality %in% c("shares_LL")){
    P <- rep(NA, el_settings$nObs)
    if(any(!el_settings$rows)) P <- apollo_insertRows(P, el_settings$rows, 1)
    return(P)
  }
  
  # ###################################################### #
  #### functionality="estimate/conditionals/raw/output" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","conditionals","raw","output")){
    P <- el_settings$probs_EL(el_settings)
    if(any(!el_settings$rows)) P <- apollo_insertRows(P, el_settings$rows, 1)
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #
  
  if(functionality=="prediction"){
    if(!apollo_inputs$silent) apollo_print('Prediction not implemented for exploded logit models.')
    return(NA)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    if(!apollo_inputs$silent) apollo_print('Gradient not implemented for exploded logit models')
    return(NA)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(el_settings$el_diagnostics(el_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(el_settings$el_diagnostics(el_settings, apollo_inputs, data =FALSE))
    return(P)
  }
    
}

