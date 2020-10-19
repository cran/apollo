#' Calculates probabilities of a Cross-nested Logit
#'
#' Calculates probabilities of a Cross-nested Logit model.
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
#'                       \item \strong{alternatives}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{avail}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item \strong{choiceVar}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{V}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{cnlNests}: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named according to nests. The lambda at the root is fixed to 1, and therefore does not need to be defined.
#'                       \item \strong{cnlStructure}: Numeric matrix. One row per nest and one column per alternative. Each element of the matrix is the alpha parameter of that (nest, alternative) pair.
#'                       \item \strong{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{componentName}: Character. Name given to model component.
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
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @export
apollo_cnl <- function(cnl_settings, functionality){
  ### Set or extract componentName
  modelType   = "CNL"
  if(is.null(cnl_settings[["componentName"]])){
    cnl_settings[["componentName"]] = ifelse(!is.null(cnl_settings[['componentName2']]),
                                             cnl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && cnl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 cnl_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, cnl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", cnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(cnl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load cnl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(cnl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded cnl_settings, restore the one received as argument
    if(is.null(tmp$V)           ) tmp$V            <- cnl_settings$V
    if(is.null(tmp$cnlNests)    ) tmp$cnlNests     <- cnl_settings$cnlNests
    if(is.null(tmp$cnlStructure)) tmp$cnlStructure <- cnl_settings$cnlStructure
    cnl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    cnl_settings <- apollo_preprocess(inputs = cnl_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print("No C++ optimisation available for CNL")
    cnl_settings$probs_CNL <- function(cnl_settings, all=FALSE){
      # Fix choiceVar if "raw" and choiceVar==NA
      cnl_settings$choiceNA = FALSE
      if(all(is.na(cnl_settings$choiceVar))){
        cnl_settings$choiceVar = cnl_settings$alternatives[1]
        cnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      cnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                               cnl_settings$V, cnl_settings$avail, SIMPLIFY=FALSE)
      # Extract chosen V or maximum V
      if(!all) VSubs <- Reduce('+', mapply("*", cnl_settings$Y, cnl_settings$V, SIMPLIFY=FALSE)) else VSubs <- do.call(pmax, cnl_settings$V)
      cnl_settings$V <- lapply(cnl_settings$V, "-", VSubs)
      rm(VSubs)
      # consider availabilities once before exponentiating (avoids issues if unavailable alternatives have attributes at 999)
      cnl_settings$V <- mapply('*', cnl_settings$V, cnl_settings$avail, SIMPLIFY = FALSE)
      # exponentiate utilities
      cnl_settings$V = lapply(X=cnl_settings$V, FUN=exp)
      # consider availabilities (it assumes eV and avail are in the same order)
      cnl_settings$V <- mapply('*', cnl_settings$V, cnl_settings$avail, SIMPLIFY = FALSE)
      # work out denominator for within nest probs
      denom_within = list()
      nests        = nrow(cnl_settings$cnlStructure)
      for(t in 1:nests){
        denom_within[[t]]=0
        for(j in 1:cnl_settings$nAlt){
          denom_within[[t]] = denom_within[[t]] + 
            (cnl_settings$cnlStructure[t,j]*cnl_settings$V[[cnl_settings$altnames[j]]])^(1/cnl_settings$cnlNests[[t]])
        }
      }
      # work out within-nest probs
      Pwithin = list()
      for(j in 1:cnl_settings$nAlt){
        Pwithin[[j]] = list()
        for(t in 1:nests){
          # includes a failsafe in denominator for empty nests, numerator will ensure ratio is 0 anyway
          Pwithin[[j]][[t]] = (cnl_settings$cnlStructure[t,j]*cnl_settings$V[[cnl_settings$altnames[j]]])^(1/cnl_settings$cnlNests[[t]])/(denom_within[[t]]+(denom_within[[t]]==0)) 
        }
      }
      # work out nest probs
      Pnest = list()
      denom_nest = 0
      for(t in 1:nests) denom_nest = denom_nest + denom_within[[t]]^cnl_settings$cnlNests[[t]]
      for(t in 1:nests) Pnest[[t]] = (denom_within[[t]]^cnl_settings$cnlNests[[t]])/denom_nest
      # work individual probs
      Palts = list()
      for(j in 1:cnl_settings$nAlt){
        Palts[[j]]=0
        for(t in 1:nests) Palts[[j]] = Palts[[j]] + Pnest[[t]]*Pwithin[[j]][[t]]
      }
      names(Palts) <- names(cnl_settings$V)
      if(!(all && cnl_settings$choiceNA)) Palts[["chosen"]] <- Reduce('+', mapply('*', cnl_settings$Y, Palts, SIMPLIFY=FALSE))
      if(!all) Palts <- Palts[["chosen"]]
      return(Palts)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && functionality %in% c("preprocess", "gradient")
    test <- test && all(sapply(cnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    cnl_settings$gradient <- FALSE
    if(test){
      cnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, cnl_settings$V)
      cnl_settings$gradient <- !is.null(cnl_settings$dV)
    }; rm(test)
    
    # Return cnl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      cnl_settings$V            <- NULL
      cnl_settings$cnlNests     <- NULL
      cnl_settings$cnlStructure <- NULL
      return(cnl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(cnl_settings$V, is.function))){
    cnl_settings$V = lapply(cnl_settings$V, function(f) if(is.function(f)) f() else f )
  }
  if(any(sapply(cnl_settings$cnlNests, is.function))){
    cnl_settings$cnlNests = lapply(cnl_settings$cnlNests, function(f) if(is.function(f)) f() else f )
  }
  if(is.function(cnl_settings$cnlStructure)) cnl_settings$cnlStructure <- cnl_settings$cnlStructure()
  cnl_settings$V <- lapply(cnl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V if neccesary
  cnl_settings$V        <- cnl_settings$V[cnl_settings$altnames]
  #cnl_settings$cnlNests <- cnl_settings$cnlNests[cnl_settings$altnames]
  if(!all(cnl_settings$rows)){
    cnl_settings$V        <- lapply(cnl_settings$V, apollo_keepRows, r=cnl_settings$rows)
    #cnl_settings$cnlNests <- lapply(cnl_settings$cnlNests, apollo_keepRows, r=cnl_settings$rows)
  } 
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(cnl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(cnl_settings, modelType, apollo_inputs)

    testL <- cnl_settings$probs_CNL(cnl_settings)
    if(any(!cnl_settings$rows)) testL <- apollo_insertRows(testL, cnl_settings$rows, 1) # insert excluded rows with value 1
    if(all(testL==0)) stop('All observations have zero probability at starting value for model component "', cnl_settings$componentName,'"')
        if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', cnl_settings$componentName,'"'))
    return(invisible(testL))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:cnl_settings$nAlt) if(length(cnl_settings$avail[[i]])==1) cnl_settings$avail[[i]] <- rep(cnl_settings$avail[[i]], cnl_settings$nObs)
    # number of available alts in each observation
    nAvAlt <- rowSums(matrix(unlist(cnl_settings$avail), ncol=cnl_settings$nAlt))
    P = 1/nAvAlt # likelihood at zero
    if(any(!cnl_settings$rows)) P <- apollo_insertRows(P, cnl_settings$rows, 1)
    return(P)
  }

  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #

  if(functionality %in% c("estimate","conditionals", "output", "components")){
    P <- cnl_settings$probs_CNL(cnl_settings, all=FALSE)
    if(any(!cnl_settings$rows)) P <- apollo_insertRows(P, cnl_settings$rows, 1) # insert excluded rows with value 1
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- cnl_settings$probs_CNL(cnl_settings, all=TRUE)
    if(any(!cnl_settings$rows)) P <- lapply(P, apollo_insertRows, r=cnl_settings$rows, val=1) # insert excluded rows with value 1
    return(P)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- utils::capture.output(apollo_diagnostics(cnl_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- utils::capture.output(apollo_diagnostics(cnl_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
}
