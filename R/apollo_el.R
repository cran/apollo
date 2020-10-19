#' Calculates Exploded Logit probabilities
#'
#' Calculates the probabilities of an Exploded Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' The function calculates the probability of a ranking as a product of Multinomial Logit models with gradually reducing availability, where scale differences can be allowed for.
#' @param el_settings List of inputs of the Exploded Logit model. It shoud contain the following.
#'                    \itemize{
#'                     \item \strong{\code{"alternatives"}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                     \item \strong{\code{"avail"}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                     \item \strong{\code{"choiceVars"}}: List of numeric vectors. Contain choices for each position of the ranking. The list must be ordered with the best choice first, second best second, etc. It will usually be a list of columns from the database. Use value -1 if a stage does not apply for a given observations (e.g. when some individuals have shorter rankings).
#'                     \item \strong{\code{"V"}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                     \item \strong{\code{"scales"}}: List of vectors. Scale factors of each Logit model. At least one element should be normalized to 1. If omitted, scale=1 for all positions is assumed.
#'                     \item \strong{\code{"rows"}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                     \item \strong{\code{"componentName"}}: Character. Name given to model component.
#'                    }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"}: Used for model estimation.
#'                        \item \code{"prediction"}: Used for model predictions.
#'                        \item \code{"validate"}: Used for validating input.
#'                        \item \code{"zero_LL"}: Used for calculating null likelihood.
#'                        \item \code{"conditionals"}: Used for calculating conditionals.
#'                        \item \code{"output"}: Used for preparing output after model estimation.
#'                        \item \code{"raw"}: Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: Not applicable (\code{NA}).
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
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
      # Loop over stages, calculating the loglikelihood for each of them
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
      # Transform loglikelihood to likelihood
      P <- exp(P)
      return(P)
    }
    
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
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(el_settings, modelType, apollo_inputs)
    
    testL <- el_settings$probs_EL(el_settings)
    testL <- apollo_insertRows(testL, el_settings$rows, 1)
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
  
  # ###################################################### #
  #### functionality="estimate/conditionals/raw/output" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","conditionals","raw","output")){
    P <- el_settings$probs_EL(el_settings)
    P <- apollo_insertRows(P, el_settings$rows, 1)
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
    P$data  <- capture.output(apollo_diagnostics(el_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(el_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
    
}

