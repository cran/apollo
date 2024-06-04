#' Calculates Ordered Logit probabilities
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
#'                     \item \strong{\code{coding}}: Numeric or character vector. Optional argument. Defines the order of the levels in \code{outcomeOrdered}. The first value is associated with the lowest level of \code{outcomeOrdered}, and the last one with the highest value. If not provided, is assumed to be \code{1:(length(tau) + 1)}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
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
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{ol_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Dependent variable overview.
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
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
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", ol_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utility by V if used
  if(!is.null(ol_settings[["utility"]])) names(ol_settings)[which(names(ol_settings)=="utility")]="V"
  
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
    ol_settings$probs_OL <- function(ol_settings, all=FALSE){
      J    <- length(ol_settings$tau) + 1
      tau1 <- Reduce("+", mapply("*", ol_settings$Y[-J], ol_settings$tau, SIMPLIFY=FALSE))
      tau0 <- Reduce("+", mapply("*", ol_settings$Y[-1], ol_settings$tau, SIMPLIFY=FALSE))
      tau1 <- apollo_setRows(tau1, ol_settings$outcomeOrdered==J,  Inf)
      tau0 <- apollo_setRows(tau0, ol_settings$outcomeOrdered==1, -Inf)
      P <- 1/(1 + exp(ol_settings$V-tau1)) - 1/(1 + exp(ol_settings$V-tau0))
      if(all){
        P2 = list()
        ol_settings$tau <- c(-Inf, ol_settings$tau, Inf)
        for(j in 1:(length(ol_settings$tau)-1)) P2[[j]] = 
            1/(1 + exp(ol_settings$V-ol_settings$tau[[j+1]])) - 1/(1 + exp(ol_settings$V-ol_settings$tau[[j]]))
        names(P2) <- ol_settings$coding
        if(!(length(ol_settings$outcomeOrdered)==1 && is.na(ol_settings$outcomeOrdered))) P2[["chosen"]] <- P
        P <- P2
      }
      return(P)
    }
    
    ol_settings$ol_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
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
    ol_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(ol_settings$V) && all(sapply(ol_settings$tau, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    ol_settings$gradient <- FALSE
    if(test){
      tmp1 <- apollo_dVdB(apollo_beta, apollo_inputs, list(V=ol_settings$V))
      tmp2 <- apollo_dVdB(apollo_beta, apollo_inputs, ol_settings$tau      )
      if(!is.null(tmp1) && !is.null(tmp2)){
        ol_settings$dV   <- lapply(tmp1, `[[`, 1) # list(dV/db1, dV/db2, ...)
        ol_settings$dTau <- tmp2 # list(b1=list(dt1/db1, dt2/db1, ...), b2=...)
      }; rm(tmp1, tmp2)
      test <- is.list(ol_settings$dV) && all(sapply(ol_settings$dV, is.function))
      test <- test && is.list(ol_settings$dTau) 
      test <- test && all(sapply(ol_settings$dTau, sapply, is.function))
      ol_settings$gradient <- test
    }; rm(test)
    
    # Construct necessary input for hessian
    test <- !is.null(ol_settings$gradient) && ol_settings$gradient && apollo_inputs$apollo_control$analyticHessian
    ol_settings$hessian <- TRUE
    if(test){
      ol_settings$ddV <- list()
      ol_settings$ddTau <- list()
      pars=names(ol_settings$dV)
      for(k1 in pars){
        ol_settings$ddV[[k1]] <- list()
        ol_settings$ddTau[[k1]] <- list()
        for(k2 in pars){
          ol_settings$ddV[[k1]][[k2]] <- Deriv::Deriv(f=ol_settings$dV[[k1]], x=k2)
          if(is.null(ol_settings$ddV[[k1]][[k2]])) ol_settings$hessian <- FALSE
          ol_settings$ddTau[[k1]][[k2]] <- list() 
          for(j in 1:length(ol_settings$tau)){  ### loop over thresholds, can do the V one outside loop (above)
            ol_settings$ddTau[[k1]][[k2]][[j]] <- Deriv::Deriv(f=ol_settings$dTau[[k1]][[j]], x=k2)
            if(is.null(ol_settings$ddTau[[k1]][[k2]][[j]])) ol_settings$hessian <- FALSE
          }
        }
      }
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
    
    if(!apollo_inputs$apollo_control$noDiagnostics) ol_settings$ol_diagnostics(ol_settings, apollo_inputs)
    testL = ol_settings$probs_OL(ol_settings)
    if(any(!ol_settings$rows)) testL <- apollo_insertRows(testL, ol_settings$rows, 1)
    if(all(testL==0)) stop("\nCALCULATION ISSUE - All observations have zero probability at starting value for model component \"",ol_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', 
                                          ol_settings$componentName,'".'), type="i")
    return(invisible(testL))
  }

  # ########################################################## #
  #### functionality="zero_LL"                              ####
  # ########################################################## #

  if(functionality=="zero_LL"){
    P <- rep(1/length(ol_settings$coding), ol_settings$nObs)
    if(any(!ol_settings$rows)) P <- apollo_insertRows(P, ol_settings$rows, 1)
    return(P)
  }

  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality %in% c("shares_LL")){
    Y       <- do.call(cbind, ol_settings$Y)
    Yshares <- colSums(Y)/nrow(Y)
    P       <- as.vector(Y%*%Yshares)
    if(any(!ol_settings$rows)) P <- apollo_insertRows(P, ol_settings$rows, 1)
    return(P)
  }
  
  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #
  
  if(functionality %in% c("estimate","conditionals","output", "components")){
    P <- ol_settings$probs_OL(ol_settings, all=FALSE)
    if(any(!ol_settings$rows)) P <- apollo_insertRows(P, ol_settings$rows, 1)
    return(P)
  } 
  
  if(functionality %in% c("prediction", "raw")){
    P <- ol_settings$probs_OL(ol_settings, all=TRUE)
    if(any(!ol_settings$rows)) P <- lapply(P, apollo_insertRows, r=ol_settings$rows, val=NA)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(ol_settings$gradient) || !ol_settings$gradient) stop("INTERNAL ISSUE - Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_ol could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_ol could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities
    Pcho <- ol_settings$probs_OL(ol_settings, all=FALSE)
    if(anyNA(Pcho)) stop("INTERNAL ISSUE - NAs in choice probabilities when calculating the gradient for component ", ol_settings$componentName)
    
    # Calculate gradient
    J  <- length(ol_settings$tau) + 1 # number of levels (alternatives)
    K  <- length(ol_settings$dV) # number of parameters
    G  <- list()
    r  <- all(ol_settings$rows) # TRUE if all rows are used (no rows excluded)
    e  <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    t0 <- Reduce("+", mapply("*", ol_settings$Y[-1], ol_settings$tau, SIMPLIFY=FALSE))
    t1 <- Reduce("+", mapply("*", ol_settings$Y[-J], ol_settings$tau, SIMPLIFY=FALSE))
    t0 <- apollo_setRows(t0, ol_settings$outcomeOrdered==1, -Inf)
    t1 <- apollo_setRows(t1, ol_settings$outcomeOrdered==J,  Inf)
    V  <- ol_settings$V
    eVT0  <- exp(V - t0)
    eVT0[is.infinite(eVT0)] <- 0
    eVT1  <- exp(V - t1)
    for(k in 1:K){
      G[[k]] <- 0
      dV    <- ol_settings$dV[[k]]; environment(dV) <- e
      dV    <- dV()
      dTau  <- ol_settings$dTau[[k]]
      dTau  <- lapply(dTau, function(f){ environment(f) <- e; return(f())})
      if(!r){ # remove rows if necessary
        dV   <- apollo_keepRows(dV, ol_settings$rows)
        dTau <- lapply(dTau, apollo_keepRows, r=ol_settings$rows)
      }
      dTau0 <- Reduce("+", mapply("*", ol_settings$Y[-1], dTau, SIMPLIFY=FALSE))
      dTau1 <- Reduce("+", mapply("*", ol_settings$Y[-J], dTau, SIMPLIFY=FALSE))
      # No need to set borders to zero as they already are.
      G[[k]] <- eVT0*(dV - dTau0)/(1 + eVT0)^2 - eVT1*(dV - dTau1)/(1 + eVT1)^2
    }; rm(dV, dTau, dTau0, dTau1, eVT0, eVT1)
    
    # Restore rows and return
    if(!all(ol_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, ol_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=ol_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############################# #
  #### functionality="hessian" ####
  # ############################# #
  if(functionality=="hessian"){

    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - ",
                                                   "apollo_mnl could not ",
                                                   "fetch apollo_beta for ",
                                                   "hessian estimation."))
    
    # Calculate probabilities
    Pcho <- ol_settings$probs_OL(ol_settings, all=FALSE)
    if(anyNA(Pcho)) stop("INTERNAL ISSUE - NAs in choice probabilities when calculating the gradient for component ", ol_settings$componentName)
    
    J  <- length(ol_settings$tau) + 1 # number of levels (alternatives)
    K  <- length(ol_settings$dV) # number of parameters
    d1L  <- list()
    d2L  <- list()
    r  <- all(ol_settings$rows) # TRUE if all rows are used (no rows excluded)
    e  <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    t0 <- Reduce("+", mapply("*", ol_settings$Y[-1], ol_settings$tau, SIMPLIFY=FALSE))
    t1 <- Reduce("+", mapply("*", ol_settings$Y[-J], ol_settings$tau, SIMPLIFY=FALSE))
    t0 <- apollo_setRows(t0, ol_settings$outcomeOrdered==1, -Inf)
    t1 <- apollo_setRows(t1, ol_settings$outcomeOrdered==J,  Inf)
    V  <- ol_settings$V
    eVT0  <- exp(V - t0)
    eVT0[is.infinite(eVT0)] <- 0
    eVT1  <- exp(V - t1)
    pars <- names(ol_settings$dV)
    
    # Calculate gradient
    for(k in 1:K){
      d1L[[k]] <- 0
      d1V    <- ol_settings$dV[[k]]; environment(d1V) <- e
      d1V    <- d1V()
      d1Tau  <- ol_settings$dTau[[k]]
      d1Tau  <- lapply(d1Tau, function(f){ environment(f) <- e; return(f())})
      if(!r){ # remove rows if necessary
        d1V   <- apollo_keepRows(d1V, ol_settings$rows)
        d1Tau <- lapply(d1Tau, apollo_keepRows, r=ol_settings$rows)
      }
      d1Tau0 <- Reduce("+", mapply("*", ol_settings$Y[-1], d1Tau, SIMPLIFY=FALSE))
      d1Tau1 <- Reduce("+", mapply("*", ol_settings$Y[-J], d1Tau, SIMPLIFY=FALSE))
      d1L[[k]] <- eVT0*(d1V - d1Tau0)/(1 + eVT0)^2 - eVT1*(d1V - d1Tau1)/(1 + eVT1)^2
    }; rm(d1V, d1Tau, d1Tau0, d1Tau1)
    
    # Restore rows and return
    if(!all(ol_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, ol_settings$rows, 1)
      d1L    <- lapply(d1L, apollo_insertRows, r=ol_settings$rows, val=0)
    }
    
    # Calculate hessian of probability of chosen alternative
    d2L <- setNames(vector(mode="list", length=K), pars)
    for(k1 in 1:K){
      d2L[[k1]] <- setNames(vector(mode="list", length=K), pars)
      d1Vk1    <- ol_settings$dV[[k1]]; environment(d1Vk1) <- e
      d1Vk1    <- d1Vk1()
      d1Tauk1  <- ol_settings$dTau[[k1]]
      d1Tauk1  <- lapply(d1Tauk1, function(f){ environment(f) <- e; return(f())})
      if(!r){ # remove rows if necessary
        d1Vk1   <- apollo_keepRows(d1Vk1, ol_settings$rows)
        d1Tauk1 <- lapply(d1Tauk1, apollo_keepRows, r=ol_settings$rows)
      }
      for(k2 in 1:k1){
        d1Vk2    <- ol_settings$dV[[k2]]; environment(d1Vk2) <- e
        d1Vk2    <- d1Vk2()
        d1Tauk2  <- ol_settings$dTau[[k2]]
        d1Tauk2  <- lapply(d1Tauk2, function(f){ environment(f) <- e; return(f())})
        d2V    <- ol_settings$ddV[[k1]][[k2]]; environment(d2V) <- e
        d2V    <- d2V()
        d2Tau  <- ol_settings$ddTau[[k1]][[k2]]
        d2Tau  <- lapply(d2Tau, function(f){ environment(f) <- e; return(f())})
        if(!r){ # remove rows if necessary
          d1Vk2   <- apollo_keepRows(d1Vk2, ol_settings$rows)
          d1Tauk2 <- lapply(d1Tauk2, apollo_keepRows, r=ol_settings$rows)
          d2V     <- apollo_keepRows(d2V, ol_settings$rows)
          d2Tau   <- lapply(d1Tauk2, apollo_keepRows, r=ol_settings$rows)
        }
        d1Tau1k1 <- Reduce("+", mapply("*", ol_settings$Y[-J], d1Tauk1, SIMPLIFY=FALSE))
        d1Tau0k1 <- Reduce("+", mapply("*", ol_settings$Y[-1], d1Tauk1, SIMPLIFY=FALSE))
        d1Tau1k2 <- Reduce("+", mapply("*", ol_settings$Y[-J], d1Tauk2, SIMPLIFY=FALSE))
        d1Tau0k2 <- Reduce("+", mapply("*", ol_settings$Y[-1], d1Tauk2, SIMPLIFY=FALSE))
        d2Tau1   <- Reduce("+", mapply("*", ol_settings$Y[-J], d2Tau, SIMPLIFY=FALSE))        
        d2Tau0   <- Reduce("+", mapply("*", ol_settings$Y[-1], d2Tau, SIMPLIFY=FALSE))
        
        deVT1k1  <- eVT1*(d1Vk1-d1Tau1k1)
        deVT0k1  <- eVT0*(d1Vk1-d1Tau0k1)
        deVT1k2  <- eVT1*(d1Vk2-d1Tau1k2)
        deVT0k2  <- eVT0*(d1Vk2-d1Tau0k2)
        d2eVT1   <- eVT1*((d1Vk1-d1Tau1k1)*(d1Vk2-d1Tau1k2)+(d2V-d2Tau1))
        d2eVT0   <- eVT0*((d1Vk1-d1Tau0k1)*(d1Vk2-d1Tau0k2)+(d2V-d2Tau0))
        
        d2L[[k1]][[k2]] <- 1/((1+eVT1)^2)*(-d2eVT1+2*(deVT1k1*deVT1k2)/(1+eVT1))+1/((1+eVT0)^2)*(+d2eVT0-2*(deVT0k1*deVT0k2)/(1+eVT0))

        # Restore rows
        if(!all(ol_settings$rows)) d2L[[k1]][[k2]] <- apollo_insertRows(d2L[[k1]][[k2]], ol_settings$rows, 0)
        # symmetry
        d2L[[k2]][[k1]] <- d2L[[k1]][[k2]]
      }
    }
    
    # Return list with everything calculated
    return(list(like = Pcho, grad=d1L, hess=d2L))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(ol_settings$ol_diagnostics(ol_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(ol_settings$ol_diagnostics(ol_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
