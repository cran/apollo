#' Calculates Multinomial Logit probabilities
#'
#' Calculates probabilities of a Multinomial Logit model.
#'
#' @param mnl_settings List of inputs of the MNL model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{\code{V}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component.
#'                     }
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
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", mnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
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
    # Do pre-processing common to most models
    mnl_settings <- apollo_preprocess(inputs = mnl_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which mnl likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation for MNL available")
    # Using R likelihood
    mnl_settings$probs_MNL=function(mnl_settings, all=FALSE, restoreRows=TRUE){
      # Fix choiceVar if "raw" and choiceVar==NA
      mnl_settings$choiceNA = FALSE
      if(all(is.na(mnl_settings$choiceVar))){
        mnl_settings$choiceVar = mnl_settings$alternatives[1]
        mnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      mnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                               mnl_settings$V, mnl_settings$avail, SIMPLIFY=FALSE)
      # if probabilities for all alternatives are requested, then P is a list
      if(all){
        if(apollo_inputs$apollo_control$subMaxV){
          ### work with subtracting the maxV
          maxV <- do.call(pmax, mnl_settings$V)
          mnl_settings$V <- lapply(mnl_settings$V, "-", maxV)
          mnl_settings$V <- lapply(X=mnl_settings$V, FUN=exp)
          mnl_settings$V <- mapply('*', mnl_settings$V, mnl_settings$avail, SIMPLIFY = FALSE)
          denom = Reduce('+',mnl_settings$V)
          P <- lapply(mnl_settings$V, "/", denom)
        } else {
          ### work with subtracting the chosenV
          chosenV <- mapply("*", mnl_settings$Y, mnl_settings$V, SIMPLIFY=FALSE)
          chosenV <- Reduce('+', chosenV)
          mnl_settings$V <- lapply(X=mnl_settings$V, "-", chosenV)
          mnl_settings$V <- lapply(X=mnl_settings$V, FUN=exp)
          # consider availabilities (it assumes V and avail are in the same order)
          mnl_settings$V <- mapply('*', mnl_settings$V, mnl_settings$avail, SIMPLIFY = FALSE)
          # calculate the denominator of the Logit probability expression
          denom = Reduce('+',mnl_settings$V)
          P <- lapply(mnl_settings$V, "/", denom)
          if(any(sapply(P, anyNA))){
            P <- lapply(P, function(p){p[is.na(p)] <- 1; return(p)} )
            sP <- Reduce("+", P)
            P <- lapply(P, "/", sP)
          }
        }
        if(!mnl_settings$choiceNA) P[["chosen"]] <- Reduce("+", mapply("*", mnl_settings$Y, P, SIMPLIFY=FALSE))
        # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
      } else { 
        ### work with subtracting the chosenV
        chosenV <- mapply("*", mnl_settings$Y, mnl_settings$V, SIMPLIFY=FALSE)
        chosenV <- Reduce('+', chosenV)
        mnl_settings$V <- lapply(X=mnl_settings$V, "-", chosenV)
        mnl_settings$V <- lapply(X=mnl_settings$V, FUN=exp)
        # consider availabilities (it assumes V and avail are in the same order)
        mnl_settings$V <- mapply('*', mnl_settings$V, mnl_settings$avail, SIMPLIFY = FALSE)
        # calculate the denominator of the Logit probability expression
        denom = Reduce('+',mnl_settings$V)
        P <- mnl_settings$chosenAvail/denom
      }
      # insert excluded rows with value 1 if only teh chosen is requested, and 0 if all
      if(any(!mnl_settings$rows) & restoreRows){
        if(is.list(P)) P <- lapply(P, apollo_insertRows, r=mnl_settings$r, val=0) else P <- apollo_insertRows(P, mnl_settings$rows, 1)
      }
      return(P)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(mnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    mnl_settings$gradient <- FALSE
    if(test){
      mnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, mnl_settings$V)
      mnl_settings$gradient <- !is.null(mnl_settings$dV)
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
  
  ### Reorder V and drop rows if neccesary
  mnl_settings$V <- mnl_settings$V[mnl_settings$altnames]
  if(!all(mnl_settings$rows)) mnl_settings$V <- lapply(mnl_settings$V, apollo_keepRows, r=mnl_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Check that alternatives are named in altcodes and V
    if(is.null(mnl_settings$altnames) || is.null(mnl_settings$altcodes) || is.null(names(mnl_settings$V))) stop("Alternatives for model component \"",mnl_settings$componentName,"\" must be named, both in 'alternatives' and 'V'.")
    
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(mnl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(mnl_settings, modelType, apollo_inputs)
    
    testL = mnl_settings$probs_MNL(mnl_settings, all=FALSE)
    if(all(testL==0)) stop("All observations have zero probability at starting value for model component \"",mnl_settings$componentName,"\"")
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
  
  # ################################################################# #
  #### functionality="estimate/prediction/conditionals/raw/output" ####
  # ################################################################# #
  
  if(functionality %in% c("estimate","conditionals", "components", "output")){
    return(mnl_settings$probs_MNL(mnl_settings, all=FALSE))
  }
  
  if(functionality %in% c("prediction","raw")){
    return(mnl_settings$probs_MNL(mnl_settings, all=TRUE))
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necesary is available
    if(is.null(mnl_settings$dV) || !all(sapply(mnl_settings$dV, is.function))) stop("Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("apollo_mnl could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("apollo_mnl could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities and derivatives of utilities for all alternatives
    P    <- mnl_settings$probs_MNL(mnl_settings, all=TRUE, restoreRows=FALSE)
    Pcho <- P[["chosen"]]
    P    <- P[-which(names(P)=="chosen")]
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    for(i in 1:length(mnl_settings$dV)) environment(mnl_settings$dV[[i]]) <- e
    dV<- lapply(mnl_settings$dV, function(dv) dv())
    if(!all(mnl_settings$rows)) for(i in 1:length(dV)) dV[[i]] <- lapply(dV[[i]], apollo_keepRows, mnl_settings$rows)
    for(i in 1:mnl_settings$nAlt) dV[[i]] <- lapply(dV[[i]], 
                                                    function(dvik){ # Make dV=0 for unavailable alternatives
                                                      test <- length(dvik)==1 && length(mnl_settings$avail[[i]])>1
                                                      if(test) dvik <- rep(dvik, mnl_settings$nObs)
                                                      dvik[!mnl_settings$avail[[i]]] <- 0
                                                      return(dvik)
                                                    })
    
    # Calculate gradient
    GA<- mapply(function(y,p) y - p, mnl_settings$Y, P, SIMPLIFY=FALSE)
    G <- list()
    for(k in 1:length(dV[[1]])){
      dVk   <- lapply(dV, function(dv) dv[[k]])
      G[[k]] <- Reduce("+", mapply("*", GA, dVk, SIMPLIFY=FALSE))
    }
    G <- lapply(G, "*", Pcho)
    
    # Restore rows and return
    if(!all(mnl_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, mnl_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=mnl_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(mnl_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(mnl_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
