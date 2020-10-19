#' Calculates the probability of an Ordered Logit model
#'
#' Calculates the probabilities of an Ordered Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' This function estimates an Ordered Logit model of the type:
#' y* = V + epsilon
#' outcomeOrdered =  1 if   -Inf < y* < tau[1]
#'      2 if tau[1] < y* < tau[2]
#'      ...
#'      maxLvl if tau[length(tau)] < y* < +Inf
#' Where epsilon is distributed standard logistic, and the values 1, 2, ..., maxLvl can be
#' replaces by coding[1], coding[2], ..., coding[maxLvl].
#' The behaviour of the function changes depending on the value of the \code{functionality} argument.
#' @param ol_settings List of settings for the OL model. It should include the following.
#'                   \itemize{
#'                     \item \code{outcomeOrdered} Numeric vector. Dependant variable. The coding of this variable is assumed to be from 1 to the maximum number of different levels. For example, if the ordered response has three possible values: "never", "sometimes" and "always", then it is assumed that outcomeOrdered contains "1" for "never", "2" for "sometimes", and 3 for "always". If another coding is used, then it should be specified using the \code{coding} argument.
#'                     \item \code{V} Numeric vector. A single explanatory variable (usually a latent variable). Must have the same number of rows as outcomeOrdered.
#'                     \item \code{tau} List of numeric vector/matrix/3-dim arrays. Thresholds. As many elements as number of different levels in the dependent variable - 1. Extreme thresholds are fixed at -inf and +inf. Mixing is allowed in thresholds.
#'                     \item \code{coding} Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{outcomeOrdered}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.
#'                     \item \code{rows} Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded. It must have length equal to the length of argument \code{outcomeOrdered}. Default value is \code{"all"}, meaning all rows are considered in the calculation.
#'                     \item \code{componentName} Character. Name given to model component.
#'                   }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"} Used for model estimation.
#'                        \item \code{"prediction"} Used for model predictions.
#'                        \item \code{"validate"} Used for validating input.
#'                        \item \code{"zero_LL"} Used for calculating null likelihood.
#'                        \item \code{"conditionals"} Used for calculating conditionals.
#'                        \item \code{"output"} Used for preparing output after model estimation.
#'                        \item \code{"raw"} Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all possible levels, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @importFrom utils capture.output
#' @export
apollo_ol  <- function(ol_settings, functionality){
  ### Set or extract componentName
  modelType = "OL"
  if(is.null(ol_settings[["componentName"]])){
    ol_settings[["componentName"]] = ifelse(!is.null(ol_settings[['componentName2']]),
                                            ol_settings[['componentName2']], modelType)
    test <- functionality=="validate" && ol_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 ol_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, ol_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", ol_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(ol_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    
    # Load ol_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(ol_settings$componentName, "_settings")]]
    # If there is no V inside the loaded ol_settings, restore the one received as argument
    if(is.null(tmp$V)  ) tmp$V <- ol_settings$V
    if(is.null(tmp$tau)) if(is.list(tmp$tau)) tmp$tau <- ol_settings$tau else tmp$tau <- as.list(ol_settings$tau)
    ol_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    ol_settings <- apollo_preprocess(inputs=ol_settings, modelType, 
                                     functionality, apollo_inputs)
    
    # Determine which mnl likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation available for OL components.")
    # Using R likelihood
    ol_settings$probs_OL <- function(ol_settings, all=FALSE, restoreRows=TRUE){
      nThr <- length(ol_settings$tau) # number of thresholds
      tau1 <- Reduce("+", mapply("*", ol_settings$Y[-(nThr+1)], ol_settings$tau, SIMPLIFY=FALSE))
      tau0 <- Reduce("+", mapply("*", ol_settings$Y[-1       ], ol_settings$tau, SIMPLIFY=FALSE))
      tau1[ol_settings$outcomeOrdered==length(ol_settings$coding)] <- Inf
      tau0[ol_settings$outcomeOrdered==1] <- -Inf
      P <- 1/(1 + exp(ol_settings$V-tau1)) - 1/(1 + exp(ol_settings$V-tau0))
      if(restoreRows && any(!ol_settings$rows)) P <- apollo_insertRows(P, ol_settings$rows, 1)
      if(all){
        P2 = list()
        ol_settings$tau <- c(-Inf, ol_settings$tau, Inf)
        for(j in 1:(length(ol_settings$tau)-1)) P2[[j]] = 
            1/(1 + exp(ol_settings$V-ol_settings$tau[[j+1]])) - 1/(1 + exp(ol_settings$V-ol_settings$tau[[j]]))
        names(P2) <- ol_settings$coding
        if(restoreRows && any(!ol_settings$rows)) P2 <- lapply(P2, apollo_insertRows, r=ol_settings$rows, val=1)
        if(!(length(ol_settings$outcomeOrdered)==1 && is.na(ol_settings$outcomeOrdered))) P2[["chosen"]] <- P
        P <- P2
      }
      return(P)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(ol_settings$V) && all(sapply(ol_settings$tau, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    ol_settings$gradient <- FALSE
    if(test){
      tmp <- list(V=ol_settings$V); tmp <- c(tmp, ol_settings$tau)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, tmp)
      if(!is.null(tmp)){
        ol_settings$dV   <- tmp$V
        ol_settings$dTau <- tmp[-which(names(tmp)=="V")]
      }
      test <- !is.null(ol_settings$dV) && is.function(ol_settings$dV)
      test <- test && !is.null(ol_settings$dTau) && is.list(ol_settings$dTau) 
      test <- test && all(sapply(ol_settings$dTau, is.function))
      ol_settings$gradient <- test
    }; rm(test)
    
    # Return ol_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      ol_settings$V   <- NULL
      ol_settings$tau <- NULL
      return(ol_settings)
    }
  }
  
  
  # ################################################## #
  #### Transform V & tau into numeric and drop rows ####
  # ################################################## #
  
  ### Evaluate V, tau
  if(is.function(ol_settings$V)) ol_settings$V <- ol_settings$V()
  if(any(sapply(ol_settings$tau, is.function))){
    ol_settings$tau <- lapply(ol_settings$tau, function(f) if(is.function(f)) f() else f)
  }
  if(is.matrix(ol_settings$V) && ncol(ol_settings$V)==1) ol_settings$V <- as.vector(ol_settings$V)
  ol_settings$tau <- lapply(ol_settings$tau, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Filter rows
  if(any(!ol_settings$rows)){
    ol_settings$V   <- apollo_keepRows(ol_settings$V, ol_settings$rows)
    ol_settings$tau <- lapply(ol_settings$tau, apollo_keepRows, r=ol_settings$rows)
  }

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(ol_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(ol_settings, modelType, apollo_inputs)
    testL = ol_settings$probs_OL(ol_settings)
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"",ol_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', 
                                          ol_settings$componentName,'".'))
    return(invisible(testL))
  }

  # ########################################################## #
  #### functionality="zero_LL"                              ####
  # ########################################################## #

  if(functionality=="zero_LL"){
    P <- rep(NA, ol_settings$nObs)
    if(any(!ol_settings$rows)) P <- apollo_insertRows(P, ol_settings$rows, 1)
    return(P)
  }

  # ################################################################## #
  #### functionality = estimate, conditional, raw, prediction & raw ####
  # ################################################################## #

  if(functionality %in% c("estimate","conditionals","output")) return(ol_settings$probs_OL(ol_settings, all=FALSE))
  
  if(functionality %in% c("prediction", "raw")) return(ol_settings$probs_OL(ol_settings, all=TRUE))
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(ol_settings$gradient) || !ol_settings$gradient) stop("Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("apollo_ol could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("apollo_ol could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities and derivatives of utilities for all alternatives
    Pcho <- ol_settings$probs_OL(ol_settings, all=FALSE, restoreRows=FALSE)
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    dV <- ol_settings$dV; environment(dV) <- e
    dV <- dV()
    for(i in 1:length(ol_settings$dTau)) environment(ol_settings$dTau[[i]]) <- e
    dTau <- lapply(ol_settings$dTau, function(dT) dT())
    if(!all(ol_settings$rows)){
      dV <- lapply(dV, apollo_keepRows, ol_settings$rows)
      for(i in 1:length(dTau)) dTau[[i]] <- lapply(dTau[[i]], apollo_keepRows, ol_settings$rows)
    }; rm(e)
    
    # Extract right thresholds
    nPar <- length(dV) # number of parameters
    nThr <- length(ol_settings$tau) # number of thresholds
    dTau1 <- list(); dTau0 <- list()
    for(i in 1:nPar){
      dTau1[[i]] <- Reduce("+", mapply("*", ol_settings$Y[-(nThr+1)], lapply(dTau, function(dT) dT[[i]]), SIMPLIFY=FALSE) )
      dTau0[[i]] <- Reduce("+", mapply("*", ol_settings$Y[-1       ], lapply(dTau, function(dT) dT[[i]]), SIMPLIFY=FALSE) )
      #dTau1[choFirst] <- Inf # No need as these will be zero be default (the derivative of Inf)
      #dTau0[choLast] <- -Inf
    }
    rm(dTau)
    tau1 <- Reduce("+", mapply("*", ol_settings$Y[-(nThr+1)], ol_settings$tau, SIMPLIFY=FALSE))
    tau0 <- Reduce("+", mapply("*", ol_settings$Y[-1       ], ol_settings$tau, SIMPLIFY=FALSE))
    tau1[ol_settings$outcomeOrdered==length(ol_settings$coding)] <- Inf
    tau0[ol_settings$outcomeOrdered==1] <- -Inf
    
    # Calculate gradient
    tmp1 <- mapply("-", dV, dTau1, SIMPLIFY=FALSE)
    tmp0 <- mapply("-", dV, dTau0, SIMPLIFY=FALSE)
    tmpA <- exp(ol_settings$V - tau1)
    tmpA <- tmpA/(1 + tmpA)^2
    tmpB <- exp(ol_settings$V - tau0)
    tmpB[tmpB==Inf] <- 0
    tmpB <- tmpB/(1 + tmpB)^2
    G    <- mapply(function(t1, t0) -tmpA*t1 + tmpB*t0, tmp1, tmp0, SIMPLIFY=FALSE)
    
    # Restore rows and return
    if(!all(ol_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, ol_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=ol_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(ol_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(ol_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
