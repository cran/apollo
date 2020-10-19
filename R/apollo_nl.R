#' Calculates probabilities of a Nested Logit
#'
#' Calculates probabilities of a Nested Logit model.
#'
#' In this implementation of the Nested Logit model, each nest must have a lambda parameter associated to it.
#' For the model to be consistent with utility maximisation, the estimated value of the Lambda parameter of all nests
#' should be between 0 and 1. Lambda parameters are inversely proportional to the correlation between the error terms of 
#' alternatives in a nest. If lambda=1, then there is no relevant correlation between the unobserved
#' utility of alternatives in that nest.
#' The tree must contain an upper nest called \code{"root"}. The lambda parameter of the root is automatically
#' set to 1 if not specified in \code{nlNests}. And while setting it to another value is possible, it is not
#' recommended.
#' @param nl_settings List of inputs of the NL model. It shoud contain the following.
#'                    \itemize{
#'                       \item \code{alternatives}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \code{avail}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item \code{choiceVar}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \code{V}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \code{nlNests}: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named with the nest name. The lambda at the root is fixed to 1 if excluded (recommended).
#'                       \item \code{nlStructure}: Named list of character vectors. As many elements as nests, it must include the "root". Each element contains the names of the nests or alternatives that belong to it. Element names must match those in \code{nlNests}.
#'                       \item \code{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \code{componentName}: Character. Name given to model component.
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
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @importFrom utils capture.output
#' @export
apollo_nl <- function(nl_settings, functionality){
  ### Set or extract componentName
  modelType   = "NL"
  if(is.null(nl_settings[["componentName"]])){
    nl_settings[["componentName"]] = ifelse(!is.null(nl_settings[['componentName2']]),
                                            nl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && nl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 nl_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, nl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", nl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(nl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load nl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(nl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded nl_settings, restore the one received as argument
    if(is.null(tmp$V)          ) tmp$V           <- nl_settings$V
    if(is.null(tmp$nlNests)    ) tmp$nlNests     <- nl_settings$nlNests
    if(is.null(tmp$nlStructure)) tmp$nlStructure <- nl_settings$nlStructure
    nl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    nl_settings <- apollo_preprocess(inputs = nl_settings, modelType, 
                                     functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print("No C++ optimisation available for NL")
    nl_settings$probs_NL <- function(nl_settings, all=FALSE){
      # Fix choiceVar if "raw" and choiceVar==NA
      nl_settings$choiceNA = FALSE
      if(all(is.na(nl_settings$choiceVar))){
        nl_settings$choiceVar = nl_settings$alternatives[1]
        nl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      nl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                              nl_settings$V, nl_settings$avail, SIMPLIFY=FALSE)
      # Extract chosen V or maximum V
      if(!all) VSubs <- Reduce('+', mapply("*", nl_settings$Y, nl_settings$V, SIMPLIFY=FALSE)) else VSubs <- do.call(pmax, nl_settings$V)
      nl_settings$V <- lapply(nl_settings$V, "-", VSubs)
      rm(VSubs)
      # Not sure what the two following lines are supposed to used for
      #combined_elements="root"
      #for(j in 1:length(nlStructure)) combined_elements=c(combined_elements,nlStructure[[j]])
      # loop over nests to create new utility elements and new availability terms
      for(k in length(nl_settings$nlStructure):1){
        nestK <- names(nl_settings$nlStructure)[k]
        nl_settings$V[[nestK]] = 0
        # calculate availability of nest
        nl_settings$avail[[nestK]] = 1*( Reduce('+', nl_settings$avail[ nl_settings$nlStructure[[k]] ]) > 0 )
        for(j in 1:length(nl_settings$nlStructure[[k]])){
          nodeJ <- nl_settings$nlStructure[[k]][j]
          nl_settings$V[[nestK]] = nl_settings$V[[nestK]] + 
            nl_settings$avail[[nodeJ]]*exp( nl_settings$V[[nodeJ]]/nl_settings$nlNests[[nestK]] )
        }
        nl_settings$V[[nestK]] = nl_settings$nlNests[[nestK]]*log(nl_settings$V[[nestK]])
      }
      # calculate log(probabilities)
      logPalts=list()
      for(j in 1:length(nl_settings$altnames)){
        logPalts[[j]]=0
        ancestorsJ <- nl_settings$ancestors[[nl_settings$altnames[[j]]]]
        for(k in 1:(length(ancestorsJ)-1)){ # loop to level just below root
          current_V = nl_settings$V[[ ancestorsJ[k] ]]
          next_V    = nl_settings$V[[ ancestorsJ[k+1] ]]
          logPalts[[j]] = logPalts[[j]] + (current_V-next_V)/nl_settings$nlNests[[ ancestorsJ[k+1] ]]
        }
      }
      Palts = lapply(X=logPalts, FUN=exp)
      names(Palts)=names(nl_settings$V)[1:length(nl_settings$altnames)]
      # consider availabilities (it assumes Palts and avail are in the same order)
      Palts <- mapply('*', Palts, nl_settings$avail[1:length(nl_settings$altnames)], SIMPLIFY = FALSE)
      Palts <- lapply(Palts, function(x) {
        x[is.na(x)] <- 0
        return(x)}) # replace all NaN by 0
      # Prepare output
      if(!(all && nl_settings$choiceNA)) Palts[["chosen"]] <- Reduce('+', mapply('*', nl_settings$Y, Palts, SIMPLIFY=FALSE))
      if(!all) Palts <- Palts[["chosen"]]
      return(Palts)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && functionality %in% c("preprocess", "gradient")
    test <- test && all(sapply(nl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    nl_settings$gradient <- FALSE
    if(test){
      nl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, nl_settings$V)
      nl_settings$gradient <- !is.null(nl_settings$dV)
    }; rm(test)
    
    # Return nl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      nl_settings$V           <- NULL
      nl_settings$nlNests     <- NULL
      nl_settings$nlStructure <- NULL
      return(nl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(nl_settings$V, is.function))){
    nl_settings$V = lapply(nl_settings$V, function(f) if(is.function(f)) f() else f )
  }
  if(any(sapply(nl_settings$nlNests, is.function))){
    nl_settings$nlNests = lapply(nl_settings$nlNests, function(f) if(is.function(f)) f() else f )
  }
  if(is.function(nl_settings$nlStructure)) nl_settings$nlStructure <- nl_settings$nlStructure()
  nl_settings$V <- lapply(nl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V if neccesary
  nl_settings$V        <- nl_settings$V[nl_settings$altnames]
  if(!all(nl_settings$rows)) nl_settings$V <- lapply(nl_settings$V, apollo_keepRows, r=nl_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(nl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(nl_settings, modelType, apollo_inputs)
    
    testL=nl_settings$probs_NL(nl_settings)
    if(any(!nl_settings$rows)) testL <- apollo_insertRows(testL, nl_settings$rows, 1) # insert excluded rows with value 1
    if(all(testL==0)) stop('All observations have zero probability at starting value for model component "', nl_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', nl_settings$componentName,'"'))
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:nl_settings$nAlt) if(length(nl_settings$avail[[i]])==1) nl_settings$avail[[i]] <- rep(nl_settings$avail[[i]], nl_settings$nObs)
    # number of available alts in each observation
    nAvAlt <- rowSums(matrix(unlist(nl_settings$avail), ncol=nl_settings$nAlt))
    P = 1/nAvAlt # likelihood at zero
    if(any(!nl_settings$rows)) P <- apollo_insertRows(P, nl_settings$rows, 1)
    return(P)
  }
  
  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #
  
  if(functionality %in% c("estimate","conditionals", "output", "components")){
    P <- nl_settings$probs_NL(nl_settings, all=FALSE)
    if(any(!nl_settings$rows)) P <- apollo_insertRows(P, nl_settings$rows, 1) # insert excluded rows with value 1
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- nl_settings$probs_NL(nl_settings, all=TRUE)
    if(any(!nl_settings$rows)) P <- lapply(P, apollo_insertRows, r=nl_settings$rows, val=1) # insert excluded rows with value 1
    return(P)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(nl_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(nl_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
