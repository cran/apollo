#' Calculates the probability of an ordered probit model
#'
#' Calculates the probabilities of an ordered probit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' This function estimates an ordered probit model of the type:
#' \deqn{ y^{*} = V + \epsilon \\
#' y = 1 if -\infty < y^{*} < \tau_1,
#'     2 if \tau_1 < y^{*} < \tau_2,
#'     ...,
#'     max(y) if \tau_{max(y)-1} < y^{*} < \infty}
#' Where \eqn{\epsilon} is distributed standard normal, and the values 1, 2, ..., \eqn{max(y)} can be
#' replaced by \code{coding[1], coding[2], ..., coding[maxLvl]}.
#' The behaviour of the function changes depending on the value of the \code{functionality} argument.
#' @param op_settings List of settings for the OP model. It should include the following.
#'                   \itemize{
#'                     \item \code{outcomeOrdered} Numeric vector. Dependant variable. The coding of this variable is assumed to be from 1 to the maximum number of different levels. For example, if the ordered response has three possible values: "never", "sometimes" and "always", then it is assumed that outcomeOrdered contains "1" for "never", "2" for "sometimes", and 3 for "always". If another coding is used, then it should be specified using the \code{coding} argument.
#'                     \item \code{V} Numeric vector/matrix/3-sim array. A single explanatory variable (usually a latent variable). Must have the same number of rows as outcomeOrdered.
#'                     \item \code{tau} List of numeric vectors/matrices/3-dim arrays. Thresholds. As many as number of different levels in the dependent variable - 1. Extreme thresholds are fixed at -inf and +inf. Mixing is allowed in thresholds. Can also be a matrix with as many rows as observations and as many columns as thresholds.
#'                     \item \code{coding} Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{outcomeOrdered}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.
#'                     \item \code{rows} Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded. It must have length equal to the length of argument \code{outcomeOrdered}. Default value is \code{"all"}, meaning all rows are considered in the calculation.
#'                     \item \code{componentName} Character. Name given to model component.
#'                   }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \strong{"estimate"} Used for model estimation.
#'                        \item \strong{"prediction"} Used for model predictions.
#'                        \item \strong{"validate"} Used for validating input.
#'                        \item \strong{"zero_LL"} Used for calculating null likelihood.n Not implemented for ordered probit.
#'                        \item \strong{"conditionals"} Used for calculating conditionals.
#'                        \item \strong{"output"} Used for preparing output after model estimation.
#'                        \item \strong{"raw"} Used for debugging.
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
#' @importFrom stats pnorm dnorm
#' @importFrom utils capture.output
#' @export
apollo_op  <- function(op_settings, functionality){
  ### Set or extract componentName
  modelType = "OP"
  if(is.null(op_settings[["componentName"]])){
    op_settings[["componentName"]] = ifelse(!is.null(op_settings[['componentName2']]),
                                            op_settings[['componentName2']], modelType)
    test <- functionality=="validate" && op_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 op_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, op_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", op_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(op_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    
    # Load op_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(op_settings$componentName, "_settings")]]
    # If there is no V inside the loaded op_settings, restore the one received as argument
    if(is.null(tmp$V)  ) tmp$V <- op_settings$V
    if(is.null(tmp$tau)) if(is.list(tmp$tau)) tmp$tau <- op_settings$tau else tmp$tau <- as.list(op_settings$tau)
    op_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    op_settings <- apollo_preprocess(inputs=op_settings, modelType, 
                                     functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation available for OP components.")
    # Using R likelihood
    op_settings$probs_OL <- function(op_settings, all=FALSE, restoreRows=TRUE){
      nThr <- length(op_settings$tau) # number of thresholds
      tau1 <- Reduce("+", mapply("*", op_settings$Y[-(nThr+1)], op_settings$tau, SIMPLIFY=FALSE))
      tau0 <- Reduce("+", mapply("*", op_settings$Y[-1       ], op_settings$tau, SIMPLIFY=FALSE))
      tau1[op_settings$outcomeOrdered==length(op_settings$coding)] <- Inf
      tau0[op_settings$outcomeOrdered==1] <- -Inf
      P <- pnorm(tau1 - op_settings$V) - pnorm(tau0 - op_settings$V)
      if(restoreRows && any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
      if(all){
        P2 = list()
        op_settings$tau <- c(-Inf, op_settings$tau, Inf)
        for(j in 1:(length(op_settings$tau)-1)) P2[[j]] = 
            pnorm(op_settings$tau[[j+1]] - op_settings$V) - pnorm(op_settings$tau[[j]] - op_settings$V)
        names(P2) <- op_settings$coding
        if(restoreRows && any(!op_settings$rows)) P2 <- lapply(P2, apollo_insertRows, r=op_settings$rows, val=1)
        if(!(length(op_settings$outcomeOrdered)==1 && is.na(op_settings$outcomeOrdered))) P2[["chosen"]] <- P
        P <- P2
      }
      return(P)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(op_settings$V) && all(sapply(op_settings$tau, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    op_settings$gradient <- FALSE
    if(test){
      tmp <- list(V=op_settings$V); tmp <- c(tmp, op_settings$tau)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, tmp)
      if(!is.null(tmp)){
        op_settings$dV   <- tmp$V
        op_settings$dTau <- tmp[-which(names(tmp)=="V")]
      }
      test <- !is.null(op_settings$dV) && is.function(op_settings$dV)
      test <- test && !is.null(op_settings$dTau) && is.list(op_settings$dTau) 
      test <- test && all(sapply(op_settings$dTau, is.function))
      op_settings$gradient <- test
    }; rm(test)
    
    # Return op_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      op_settings$V   <- NULL
      op_settings$tau <- NULL
      return(op_settings)
    }
  }
  
  
  # ################################################## #
  #### Transform V & tau into numeric and drop rows ####
  # ################################################## #
  
  ### Evaluate V, tau
  if(is.function(op_settings$V)) op_settings$V <- op_settings$V()
  if(any(sapply(op_settings$tau, is.function))){
    op_settings$tau <- lapply(op_settings$tau, function(f) if(is.function(f)) f() else f)
  }
  if(is.matrix(op_settings$V) && ncol(op_settings$V)==1) op_settings$V <- as.vector(op_settings$V)
  op_settings$tau <- lapply(op_settings$tau, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Filter rows
  if(any(!op_settings$rows)){
    op_settings$V   <- apollo_keepRows(op_settings$V, op_settings$rows)
    op_settings$tau <- lapply(op_settings$tau, apollo_keepRows, r=op_settings$rows)
  }
  

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(op_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(op_settings, modelType, apollo_inputs)
    testL = op_settings$probs_OL(op_settings)
    if(all(testL==0)) stop('All observations have zero probability at starting value for model component "',op_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "',
                                          op_settings$componentName,'".'))
    return(invisible(testL))
  }

  # ########################################################## #
  #### functionality="zero_LL"                              ####
  # ########################################################## #

  if(functionality=="zero_LL"){
    P <- rep(NA, op_settings$nObs)
    if(any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
    return(P)
  }

  # ############################################################# #
  #### functionality = estimate, conditional, raw & prediction ####
  # ############################################################# #

  if(functionality %in% c("estimate","conditionals", "output")) return(op_settings$probs_OL(op_settings, all=FALSE))
  
  if(functionality %in% c("prediction", "raw")      ) return(op_settings$probs_OL(op_settings, all=TRUE))
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(op_settings$gradient) || !op_settings$gradient) stop("Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("apollo_ol could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("apollo_ol could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities and derivatives of utilities for all alternatives
    Pcho <- op_settings$probs_OL(op_settings, all=FALSE, restoreRows=FALSE)
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    dV <- op_settings$dV; environment(dV) <- e
    dV <- dV()
    for(i in 1:length(op_settings$dTau)) environment(op_settings$dTau[[i]]) <- e
    dTau <- lapply(op_settings$dTau, function(dT) dT())
    if(!all(op_settings$rows)){
      dV <- lapply(dV, apollo_keepRows, op_settings$rows)
      for(i in 1:length(dTau)) dTau[[i]] <- lapply(dTau[[i]], apollo_keepRows, op_settings$rows)
    }; rm(e)
    
    # Extract right thresholds
    nPar <- length(dV) # number of parameters
    nThr <- length(op_settings$tau) # number of thresholds
    dTau1 <- list(); dTau0 <- list()
    for(i in 1:nPar){
      dTau1[[i]] <- Reduce("+", mapply("*", op_settings$Y[-(nThr+1)], lapply(dTau, function(dT) dT[[i]]), SIMPLIFY=FALSE) )
      dTau0[[i]] <- Reduce("+", mapply("*", op_settings$Y[-1       ], lapply(dTau, function(dT) dT[[i]]), SIMPLIFY=FALSE) )
    }
    rm(dTau)
    tau1 <- Reduce("+", mapply("*", op_settings$Y[-(nThr+1)], op_settings$tau, SIMPLIFY=FALSE))
    tau0 <- Reduce("+", mapply("*", op_settings$Y[-1       ], op_settings$tau, SIMPLIFY=FALSE))
    tau1[op_settings$outcomeOrdered==length(op_settings$coding)] <- Inf
    tau0[op_settings$outcomeOrdered==1] <- -Inf
    
    # Calculate gradient
    tmpA1 <- dnorm(tau1 - op_settings$V)
    tmpB1 <- dnorm(tau0 - op_settings$V)
    tmpA2 <- mapply("-", dTau1, dV, SIMPLIFY=FALSE)
    tmpB2 <- mapply("-", dTau0, dV, SIMPLIFY=FALSE)
    G     <- mapply(function(ta2, tb2) tmpA1*ta2 - tmpB1*tb2, 
                    tmpA2, tmpB2, SIMPLIFY=FALSE)
    
    # Restore rows and return
    if(!all(op_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, op_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=op_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(op_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(op_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
