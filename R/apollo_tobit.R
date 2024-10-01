#' Calculates density for a Tobit model (censored Normal)
#'
#' Calculates density for a censored Normal distribution at a specific value with a specified mean and standard deviation and user provided bounds, and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' This function calculates the probability of the linear model outcomeTobit = mu + xTobit + epsilon, where epsilon is a random error distributed Normal(0,sigma), but with optional lower and upper bounds imposed by the user (outside of which the density would be 0).
#' @param tobit_settings List of arguments to the functions. It must contain the following.
#'                               \itemize{
#'                                 \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                                 \item \strong{\code{lowerLimit}}: Numeric scalar. Lower bound beyond which the density is 0. If not provided by the user, this will be set to -Inf.
#'                                 \item \strong{\code{mu}}: Numeric scalar. Intercept of the linear model.
#'                                 \item \strong{\code{outcomeTobit}}: Numeric vector. Dependent variable.
#'                                 \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                 \item \strong{\code{sigma}}: Numeric scalar. Variance of error component of linear model to be estimated.
#'                                 \item \strong{\code{upperLimit}}: Numeric scalar. Upper bound beyond which the density is 0. If not provided by the user, this will be set to +Inf.
#'                                 \item \strong{\code{xTobit}}: Numeric vector. Single explanatory variable.
#'                               }
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
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the likelihood for each observation.
#'           \item \strong{\code{"gradient"}}: List containing the likelihood and gradient of the model component.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: Predicted value at the observation level.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{tobit_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"report"}}: Dependent variable overview.
#'           \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'         }
#' @importFrom utils capture.output
#' @export
apollo_tobit <- function(tobit_settings, functionality){
  ### Set or extract componentName
  modelType = "Tobit"
  if(is.null(tobit_settings[["componentName"]])){
    tobit_settings[["componentName"]] = ifelse(!is.null(tobit_settings[['componentName2']]),
                                                       tobit_settings[['componentName2']], modelType)
    test <- functionality=="validate" && tobit_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 tobit_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, tobit_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", tobit_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  if( !is.null(apollo_inputs[[paste0(tobit_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(tobit_settings$componentName, "_settings")]]
    # Restore variable elements if they do not exist
    if(is.null(tmp$xTobit)) tmp$xTobit <- tobit_settings$xTobit
    if(is.null(tmp$mu     )) tmp$mu      <- tobit_settings$mu     
    if(is.null(tmp$sigma  )) tmp$sigma   <- tobit_settings$sigma  
    tobit_settings <- tmp
    rm(tmp)
  } else {
    # Do pre-processing common to most models
    tobit_settings = apollo_preprocess(tobit_settings, modelType, 
                                               functionality, apollo_inputs)
    
    
    
    
    # probability function
    tobit_settings$probs_Tobit=function(tobit_settings){
      P=(tobit_settings$below*(stats::pnorm(tobit_settings$lowerLimit-tobit_settings$xTobit,
                      mean=tobit_settings$mu,sd=tobit_settings$sigma))
       +tobit_settings$within*(stats::dnorm(tobit_settings$outcomeTobit - tobit_settings$xTobit,
                      mean=tobit_settings$mu,sd=tobit_settings$sigma))
       +tobit_settings$above*(1-stats::pnorm(tobit_settings$upperLimit-tobit_settings$xTobit,
                        mean=tobit_settings$mu,sd=tobit_settings$sigma))
      )
      return(P)
    }
    
    # diagnostics function
    tobit_settings$tobit_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      if(!apollo_inputs$silent & data){
        apollo_print('\n')
        apollo_print(paste0('Summary statistics for ', toupper(inputs$modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        tmp <- t(as.matrix(summary(inputs$outcomeTobit))); rownames(tmp) <- ""
        print(round(tmp,4))
        rm(tmp)
      }
      return(invisible(TRUE))
    }
    # Store model type
    tobit_settings$modelType <- modelType
    # Construct necessary input for analytic gradient
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    vParNames <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_inputs$apollo_fixed)]
    K <- length(vParNames)
    test <- !is.null(apollo_beta) && apollo_inputs$apollo_control$analyticGrad 
    test <- test && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(tobit_settings$xTobit) 
    test <- test && is.function(tobit_settings$mu)
    test <- test && is.function(tobit_settings$sigma)
    tobit_settings$gradient <- FALSE
    if(test){
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(tobit_settings$xTobit))
      tobit_settings$dX <- lapply(tmp, `[[`, 1)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(tobit_settings$mu     ))
      tobit_settings$dM <- lapply(tmp, `[[`, 1)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(tobit_settings$sigma  ))
      tobit_settings$dS <- lapply(tmp, `[[`, 1)
      tmp <-        all( sapply(tobit_settings$dX    , is.function) )
      tmp <- tmp && all( sapply(tobit_settings$dMu   , is.function) )
      tmp <- tmp && all( sapply(tobit_settings$dSigma, is.function) )
      tobit_settings$gradient <- tmp
    }; rm(test)
    # Construct necessary input for analytic hessian
    #test <- !is.null(apollo_beta) && apollo_inputs$apollo_control$analyticHessian
    #test <- test && (functionality %in% c("preprocess", "hessian"))
    #tobit_settings$hessian <- test
    #if(test){
    #  ddX <- setNames(vector(mode="list", length=K), vParNames)
    #  ddM <- setNames(vector(mode="list", length=K), vParNames)
    #  ddS <- setNames(vector(mode="list", length=K), vParNames)
    #  for(k1 in vParNames){
    #    ddX[[k1]] <- setNames(vector(mode="list", length=K), vParNames)
    #    ddM[[k1]] <- setNames(vector(mode="list", length=K), vParNames)
    #    ddS[[k1]] <- setNames(vector(mode="list", length=K), vParNames)
    #    for(k2 in vParNames){
    #      ddX[[k1]][[k2]] <- Deriv::Deriv(f=tobit_settings$dX[[k1]], x=k2)
    #      ddM[[k1]][[k2]] <- Deriv::Deriv(f=tobit_settings$dM[[k1]], x=k2)
    #      ddS[[k1]][[k2]] <- Deriv::Deriv(f=tobit_settings$dS[[k1]], x=k2)
    #    }
    #    tobit_settings$hessian <- !any(sapply(ddX[[k1]], is.null))
    #    tobit_settings$hessian <- !any(sapply(ddM[[k1]], is.null))
    #    tobit_settings$hessian <- !any(sapply(ddS[[k1]], is.null))
    #  }
    #  tobit_settings$ddX <- ddX
    #  tobit_settings$ddM <- ddM
    #  tobit_settings$ddS <- ddS
    #}
    
    # Return tobit_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      tobit_settings$xTobit <- NULL
      tobit_settings$mu      <- NULL
      tobit_settings$sigma   <- NULL
      return(tobit_settings)
    }
  }
  
  
  # ################################################## #
  #### Evaluate xTobit, mu & sigma, and remove rows ####
  # ################################################## #
  
  ### Evaluate xTobit, mu & sigma
  if(is.function(tobit_settings$xTobit)) tobit_settings$xTobit <- tobit_settings$xTobit()
  if(is.function(tobit_settings$mu     )) tobit_settings$mu      <- tobit_settings$mu()
  if(is.function(tobit_settings$sigma  )) tobit_settings$sigma   <- tobit_settings$sigma()
  
  ### Drop excluded rows
  if(any(!tobit_settings$rows)){
    tobit_settings$outcomeTobit <- apollo_keepRows(tobit_settings$outcomeTobit, tobit_settings$rows)
    tobit_settings$xTobit       <- apollo_keepRows(tobit_settings$xTobit      , tobit_settings$rows)
    tobit_settings$mu           <- apollo_keepRows(tobit_settings$mu          , tobit_settings$rows)
    tobit_settings$sigma        <- apollo_keepRows(tobit_settings$sigma       , tobit_settings$rows)
    tobit_settings$below        <- apollo_keepRows(tobit_settings$below       , tobit_settings$rows)
    tobit_settings$above        <- apollo_keepRows(tobit_settings$above       , tobit_settings$rows)
    tobit_settings$within       <- apollo_keepRows(tobit_settings$within      , tobit_settings$rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  
  if(functionality=="validate"){

    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
                                error=function(e) return(list(noValidation=FALSE, noDiagnostics=FALSE)) )

    if(!apollo_control$noValidation) apollo_validate(tobit_settings, modelType,
                                                     functionality, apollo_inputs)
    if(!apollo_control$noDiagnostics) tobit_settings$tobit_diagnostics(tobit_settings, apollo_inputs)
    testL = tobit_settings$probs_Tobit(tobit_settings)
    if(any(!tobit_settings$rows)) testL <- apollo_insertRows(testL, tobit_settings$rows, 1)
    if(all(testL==0)) stop("\nCALCULATION ISSUE - All observations have zero probability at starting value for model component \"", tobit_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("\nSome observations have zero probability at starting value for model component \"", tobit_settings$componentName,"\""), type="i")
    return(invisible(testL))
  }

  # ####################################### #
  #### functionality="zero_LL/shares_LL" ####
  # ####################################### #
  
  if(functionality %in% c("zero_LL","shares_LL")){
    P <- rep(NA, tobit_settings$nObs)
    if(any(!tobit_settings$rows)) P <- apollo_insertRows(P, tobit_settings$rows, 1)
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction"){
    ans <- tobit_settings$mu + tobit_settings$xTobit
    ans <- ((ans<=tobit_settings$lowerLimit)*tobit_settings$lowerLimit
            +((ans>tobit_settings$lowerLimit)*(ans<tobit_settings$upperLimit))*ans
            +(ans>=tobit_settings$upperLimit)*tobit_settings$upperLimit)
    if(any(!tobit_settings$rows)) ans <- apollo_insertRows(ans, tobit_settings$rows, NA)
    return(ans)
  } 

  # ################################################################# #
  #### functionality="estimate/conditionals/raw/output/components" ####
  # ################################################################# #

  if(functionality %in% c("estimate", "conditionals", "raw", "output", "components")){
    ans <- tobit_settings$probs_Tobit(tobit_settings)
    if(any(!tobit_settings$rows)){
      if(functionality=="raw"){
        ans <- apollo_insertRows(ans, tobit_settings$rows, NA)
      } else {
        ans <- apollo_insertRows(ans, tobit_settings$rows, 1)
      } 
    }
    return(ans)
  }
  
  # ################################ #
  #### functionality="gradient"   ####
  # ################################ #
  
  if(functionality=="gradient"){
    if(!tobit_settings$gradient) stop("INTERNAL ISSUE - Analytical gradient could not be calculated for ", 
                                              tobit_settings$componentName, 
                                              ". Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_tobit could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_tobit could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate likelihood
    L <- tobit_settings$probs_Tobit(tobit_settings)
    
    
    #stats::dnorm(tobit_settings$outcomeTobit - tobit_settings$xTobit,
    #             tobit_settings$mu, tobit_settings$sigma)
    
    # Calculate dL/db
    G <- list()
    K <- length(tobit_settings$dX) # number of parameters
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    eHatSig <- tobit_settings$outcomeTobit - tobit_settings$xTobit - tobit_settings$mu
    eHatSig <- eHatSig/tobit_settings$sigma
    r <- all(tobit_settings$rows)
    z_below <- (tobit_settings$lowerLimit - tobit_settings$xTobit - tobit_settings$mu)/tobit_settings$sigma
    z_above <- (tobit_settings$upperLimit - tobit_settings$xTobit - tobit_settings$mu)/tobit_settings$sigma
    phi_z_below <- stats::dnorm(z_below)
    phi_z_above <- stats::dnorm(z_above)
    for(k in 1:K){
      # Evaluate partial derivatives
      dX <- tobit_settings$dX[[k]]; environment(dX) <- e; dX <- dX()
      dM <- tobit_settings$dM[[k]]; environment(dM) <- e; dM <- dM()
      dS <- tobit_settings$dS[[k]]; environment(dS) <- e; dS <- dS()
      if(!r){ # remove excluded rows
        dX <- apollo_keepRows(dX, tobit_settings$rows)
        dM <- apollo_keepRows(dM, tobit_settings$rows)
        dS <- apollo_keepRows(dS, tobit_settings$rows)
      }
      G[[k]] <- tobit_settings$within * L/tobit_settings$sigma*( eHatSig*(dX + dM) + dS*(eHatSig^2 - 1) )
      if(any(tobit_settings$below!=0)) G[[k]] <- G[[k]] + tobit_settings$below * phi_z_below * ((-dX-dM)/tobit_settings$sigma-((tobit_settings$lowerLimit - tobit_settings$xTobit - tobit_settings$mu)*dS)/(tobit_settings$sigma^2))
      if(any(tobit_settings$above!=0)) G[[k]] <- G[[k]] - tobit_settings$above * phi_z_above * ((-dX-dM)/tobit_settings$sigma-((tobit_settings$upperLimit - tobit_settings$xTobit - tobit_settings$mu)*dS)/(tobit_settings$sigma^2))
      
       }; rm(dX, dM, dS)
    
    # Restore rows
    if(!r){
      L <- apollo_insertRows(L, tobit_settings$rows, 1)
      G <- lapply(G, apollo_insertRows, r=tobit_settings$rows, val=0)
    }
    return(list(like=L, grad=G))
  }
  
  
  # ############ #
  #### Hessian ####
  # ############ #
  # if(functionality=="hessian"){
  #   # Checks
  #   if(!tobit_settings$hessian) stop("INTERNAL ISSUE - Analytical hessian could not be calculated for ", 
  #                                             tobit_settings$componentName, 
  #                                             ". Please set apollo_control$analyticHessian=FALSE.")
  #   apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
  #                           error=function(e) stop("INTERNAL ISSUE - apollo_tobit could not fetch apollo_beta for gradient estimation."))
  #   if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_tobit could not fetch apollo_inputs$database for gradient estimation.")
  #   
  #   
  #   # Calculate likelihood
  #   L <- stats::dnorm(tobit_settings$outcomeTobit - tobit_settings$xTobit,
  #                     tobit_settings$mu, tobit_settings$sigma)
  #   
  #   # Calculate gradient and Hessian
  #   K <- length(tobit_settings$dX) # number of parameters
  #   parNames <- names(tobit_settings$dX)
  #   e     <- list2env(c(as.list(apollo_beta), apollo_inputs$database, 
  #                       list(apollo_inputs=apollo_inputs)), hash=TRUE)
  #   s     <- tobit_settings$sigma
  #   eHat  <- tobit_settings$outcomeTobit - 
  #     tobit_settings$xTobit - tobit_settings$mu
  #   eHats <- eHat/s
  #   r  <- all(tobit_settings$rows)
  #   G  <- setNames(vector(mode="list", length=K), parNames)
  #   H  <- setNames(vector(mode="list", length=K), parNames)
  #   dX <- setNames(vector(mode="list", length=K), parNames)
  #   dM <- setNames(vector(mode="list", length=K), parNames)
  #   dS <- setNames(vector(mode="list", length=K), parNames)
  #   for(i in 1:K){
  #     # Evaluate partial derivatives
  #     dX[[i]] <- tobit_settings$dX[[i]]; environment(dX[[i]]) <- e; dX[[i]] <- dX[[i]]()
  #     dM[[i]] <- tobit_settings$dM[[i]]; environment(dM[[i]]) <- e; dM[[i]] <- dM[[i]]()
  #     dS[[i]] <- tobit_settings$dS[[i]]; environment(dS[[i]]) <- e; dS[[i]] <- dS[[i]]()
  #     if(!r){
  #       dX[[i]] <- apollo_keepRows(dX[[i]], r=tobit_settings$rows)
  #       dM[[i]] <- apollo_keepRows(dM[[i]], r=tobit_settings$rows)
  #       dS[[i]] <- apollo_keepRows(dS[[i]], r=tobit_settings$rows)
  #     }
  #     G[[i]] <- L/s*( eHats*(dX[[i]] + dM[[i]]) + dS[[i]]*(eHats^2 - 1) )
  #     H[[i]] <- setNames(vector(mode="list", length=K), parNames)
  #     for(j in 1:i){
  #       ddX <- tobit_settings$ddX[[i]][[j]]; environment(ddX) <- e; ddX <- ddX()
  #       ddM <- tobit_settings$ddM[[i]][[j]]; environment(ddM) <- e; ddM <- ddM()
  #       ddS <- tobit_settings$ddS[[i]][[j]]; environment(ddS) <- e; ddS <- ddS()
  #       if(!r){ # remove excluded rows
  #         ddX <- apollo_keepRows(ddX, tobit_settings$rows)
  #         ddM <- apollo_keepRows(ddM, tobit_settings$rows)
  #         ddS <- apollo_keepRows(ddS, tobit_settings$rows)
  #       }
  #       H[[i]][[j]] <- (G[[j]]/L - dS[[j]]/s)*G[[i]] + 
  #         (L/s) * (eHats*(ddX + ddM) + ddS*(eHats^2 - 1) - 
  #                    ((dX[[j]] + dM[[j]])/s + eHats/s*dS[[j]])*(dX[[i]] + dM[[i]] + 2*dS[[i]]*eHats) )
  #     }
  #   }; rm(dX, dM, dS)
  #   
  #   # Restore rows and return
  #   if(!r){
  #     L <- apollo_insertRows(L, tobit_settings$rows, 1)
  #     G <- lapply(G, apollo_insertRows, r=tobit_settings$rows, val=0)
  #     for(i in 1:length(H)) H[[i]] <- lapply(H[[i]], apollo_insertRows, r=tobit_settings$rows, val=0)
  #   }
  #   return(list(like=L, grad=G, hess=H))
  # }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(tobit_settings$tobit_diagnostics(tobit_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(tobit_settings$tobit_diagnostics(tobit_settings, apollo_inputs, data =FALSE))
    return(P)
  }

}
