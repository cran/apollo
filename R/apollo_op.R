#' Calculates Ordered Probit probabilities
#'
#' Calculates the probabilities of an Ordered Probit model and can also perform other operations based on the value of the \code{functionality} argument.
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
#'                     \item \strong{\code{coding}}: Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{outcomeOrdered}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.
#'                     \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                     \item \strong{\code{outcomeOrdered}}: Numeric vector. Dependent variable. The coding of this variable is assumed to be from 1 to the maximum number of different levels. For example, if the ordered response has three possible values: "never", "sometimes" and "always", then it is assumed that outcomeOrdered contains "1" for "never", "2" for "sometimes", and 3 for "always". If another coding is used, then it should be specified using the \code{coding} argument.
#'                     \item \strong{\code{rows}}: Boolean vector. TRUE if a row must be considered in the calculations, FALSE if it must be excluded. It must have length equal to the length of argument \code{outcomeOrdered}. Default value is \code{"all"}, meaning all rows are considered in the calculation.
#'                     \item \strong{\code{tau}}: List of numeric vectors/matrices/3-dim arrays. Thresholds. As many as number of different levels in the dependent variable - 1. Extreme thresholds are fixed at -inf and +inf. Mixing is allowed in thresholds. Can also be a matrix with as many rows as observations and as many columns as thresholds.
#'                     \item \strong{\code{utilities}}: Numeric vector/matrix/3-sim array. A single explanatory variable (usually a latent variable). Must have the same number of rows as outcomeOrdered.
#'                   }
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
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all possible levels, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{op_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Dependent variable overview.
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
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
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", op_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utility by V if used
  if(!is.null(op_settings[["utility"]])) names(op_settings)[which(names(op_settings)=="utility")]="V"
  
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
    op_settings$probs_OP <- function(op_settings, all=FALSE){#}, restoreRows=TRUE){
      nThr <- length(op_settings$tau) # number of thresholds
      tau1 <- Reduce("+", mapply("*", op_settings$Y[-(nThr+1)], op_settings$tau, SIMPLIFY=FALSE))
      tau0 <- Reduce("+", mapply("*", op_settings$Y[-1       ], op_settings$tau, SIMPLIFY=FALSE))
      tau1[op_settings$outcomeOrdered==length(op_settings$coding)] <- Inf
      tau0[op_settings$outcomeOrdered==1] <- -Inf
      P <- pnorm(tau1 - op_settings$V) - pnorm(tau0 - op_settings$V)
      #if(restoreRows && any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
      if(all){
        P2 = list()
        op_settings$tau <- c(-Inf, op_settings$tau, Inf)
        for(j in 1:(length(op_settings$tau)-1)) P2[[j]] = 
            pnorm(op_settings$tau[[j+1]] - op_settings$V) - pnorm(op_settings$tau[[j]] - op_settings$V)
        names(P2) <- op_settings$coding
        #if(restoreRows && any(!op_settings$rows)) P2 <- lapply(P2, apollo_insertRows, r=op_settings$rows, val=1)
        if(!(length(op_settings$outcomeOrdered)==1 && is.na(op_settings$outcomeOrdered))) P2[["chosen"]] <- P
        P <- P2
      }
      return(P)
    }
    
    op_settings$op_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      choicematrix <- t(as.matrix(table(inputs$outcomeOrdered)))
      choicematrix <- rbind(choicematrix, choicematrix[1,]/inputs$nObs*100)
      rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")
      if(!apollo_inputs$silent & data){
        apollo_print("\n")
        apollo_print(paste0('Overview of choices for ', toupper(inputs$modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        print(round(choicematrix,2))
      }
      return(invisible(TRUE))
    }
    
    # Store model type
    op_settings$modelType <- modelType
    ### added 26 April, now in preprocess
    ### op_settings$nAlt = length(op_settings$tau) + 1
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(op_settings$V) && all(sapply(op_settings$tau, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    op_settings$gradient <- FALSE
    if(test){
      tmp1 <- apollo_dVdB(apollo_beta, apollo_inputs, list(V=op_settings$V))
      tmp2 <- apollo_dVdB(apollo_beta, apollo_inputs, op_settings$tau      )
      if(!is.null(tmp1) && !is.null(tmp2)){
        op_settings$dV   <- lapply(tmp1, `[[`, 1) # list(dV/db1, dV/db2, ...)
        op_settings$dTau <- tmp2 # list(b1=list(dt1/db1, dt2/db1, ...), b2=...)
      }; rm(tmp1, tmp2)
      test <- is.list(op_settings$dV) && sapply(op_settings$dV, is.function)
      test <- test && is.list(op_settings$dTau) 
      test <- test && all(sapply(op_settings$dTau, sapply, is.function))
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
    
    if(!apollo_inputs$apollo_control$noDiagnostics) op_settings$op_diagnostics(op_settings, apollo_inputs)
    testL = op_settings$probs_OP(op_settings)
    if(any(!op_settings$rows)) testL <- apollo_insertRows(testL, op_settings$rows, 1)
    if(all(testL==0)) stop('CALCULATION ISSUE - All observations have zero probability at starting value for model component "',op_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "',
                                          op_settings$componentName,'".'), type="i")
    return(invisible(testL))
  }

  # ########################################################## #
  #### functionality="zero_LL"                              ####
  # ########################################################## #
  
  if(functionality=="zero_LL"){
    P <- rep(1/length(op_settings$coding), op_settings$nObs)
    if(any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
    return(P)
  }
  
  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality %in% c("shares_LL")){
    Y       <- do.call(cbind, op_settings$Y)
    Yshares <- colSums(Y)/nrow(Y)
    P       <- as.vector(Y%*%Yshares)
    if(any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
    return(P)
  }
  
  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #

  if(functionality %in% c("estimate","conditionals","output", "components")){
    P <- op_settings$probs_OP(op_settings, all=FALSE)
    if(any(!op_settings$rows)) P <- apollo_insertRows(P, op_settings$rows, 1)
    return(P)
  } 
  
  if(functionality %in% c("prediction", "raw")){
    P <- op_settings$probs_OP(op_settings, all=TRUE)
    if(any(!op_settings$rows)) P <- lapply(P, apollo_insertRows, r=op_settings$rows, val=NA)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(op_settings$gradient) || !op_settings$gradient) stop("INTERNAL ISSUE - Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_ol could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_ol could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities
    Pcho <- op_settings$probs_OP(op_settings, all=FALSE)
    if(anyNA(Pcho)) stop("INTERNAL ISSUE - NAs in choice probabilities when calculating the gradient for component ", op_settings$componentName)
    
    # Calculate gradient
    J  <- length(op_settings$tau) + 1 # number of levels (alternatives)
    K  <- length(op_settings$dV) # number of parameters
    G  <- list()
    r  <- all(op_settings$rows) # TRUE if all rows are used (no rows excluded)
    e  <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    t0 <- Reduce("+", mapply("*", op_settings$Y[-1], op_settings$tau, SIMPLIFY=FALSE))
    t1 <- Reduce("+", mapply("*", op_settings$Y[-J], op_settings$tau, SIMPLIFY=FALSE))
    t0 <- apollo_setRows(t0, op_settings$outcomeOrdered==1, -Inf)
    t1 <- apollo_setRows(t1, op_settings$outcomeOrdered==J,  Inf)
    V  <- op_settings$V
    for(k in 1:K){
      G[[k]] <- 0
      dV    <- op_settings$dV[[k]]; environment(dV) <- e
      dV    <- dV()
      dTau  <- op_settings$dTau[[k]]
      dTau  <- lapply(dTau, function(f){ environment(f) <- e; return(f())})
      if(!r){ # remove rows if necessary
        dV   <- apollo_keepRows(dV, op_settings$rows)
        dTau <- lapply(dTau, apollo_keepRows, r=op_settings$rows)
      }
      dTau0 <- Reduce("+", mapply("*", op_settings$Y[-1], dTau, SIMPLIFY=FALSE))
      dTau1 <- Reduce("+", mapply("*", op_settings$Y[-J], dTau, SIMPLIFY=FALSE))
      # No need to set borders to zero as they already are.
      eVT0  <- exp(V - t0)
      eVT0[is.infinite(eVT0)] <- 0
      eVT1  <- exp(V - t1)
      G[[k]] <- eVT0*(dV - dTau0)/(1 + eVT0)^2 - eVT1*(dV - dTau1)/(1 + eVT1)^2
    }; rm(dV, dTau, dTau0, dTau1, eVT0, eVT1)
    
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
    P$data  <- capture.output(op_settings$op_diagnostics(op_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(op_settings$op_diagnostics(op_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
