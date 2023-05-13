#' Calculates Fractional Multinomial Logit probabilities
#'
#' Calculates the probabilities of a Fractional Multinomial Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' @param fmnl_settings List of inputs of the FMNL model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'utilities'.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{choiceShares}}: Named list of numeric vectors. Share allocated to each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                     }
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
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{mfnl_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Overview of dependent variable
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @export
#' @importFrom utils capture.output
apollo_fmnl <- function(fmnl_settings, functionality){
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE, silent=FALSE, analyticGrad=TRUE)) ))
  
  ### Set or extract componentName
  modelType   = "FMNL"
  if(is.null(fmnl_settings[["componentName"]])){
    fmnl_settings[["componentName"]] = ifelse(!is.null(fmnl_settings[['componentName2']]),
                                             fmnl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && fmnl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', fmnl_settings[["componentName"]], '" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, fmnl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", fmnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  ### replace utilities by V if used
  if(!is.null(fmnl_settings[["utilities"]])) names(fmnl_settings)[which(names(fmnl_settings)=="utilities")]="V"
  
  ### turn character vector of alternatives into named vector (to make compatible with MNL code)
  if(is.character(fmnl_settings[["alternatives"]])) fmnl_settings[["alternatives"]]=setNames(1:length(fmnl_settings[["alternatives"]]), fmnl_settings[["alternatives"]])
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  if( !is.null(apollo_inputs[[paste0(fmnl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load fmnl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(fmnl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded fmnl_settings, restore the one received as argument
    if(is.null(tmp$V)) tmp$V <- fmnl_settings$V
    fmnl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    fmnl_settings <- apollo_preprocess(inputs = fmnl_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which fmnl likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation for fmnl available")
    # Using R likelihood
    fmnl_settings$probs_fmnl=function(fmnl_settings, all=FALSE){#}, restoreRows=TRUE){
      # Fix choiceVar if "raw" and choiceVar==NA
      fmnl_settings$choiceNA = FALSE
      if(all(is.na(unlist(fmnl_settings$choiceShares)))){
        for(s in 1:length(fmnl_settings$choiceShares)) fmnl_settings$choiceShares[[s]] = rep(1/length(fmnl_settings$choiceShares),length(fmnl_settings$alternatives[1]))
        fmnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      fmnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                               fmnl_settings$V, fmnl_settings$avail, SIMPLIFY=FALSE)
      # if probabilities for all alternatives are requested, then P is a list
      if(all){
        if(apollo_inputs$apollo_control$subMaxV){
          ### work with subtracting the maxV
          maxV <- do.call(pmax, fmnl_settings$V)
          fmnl_settings$V <- lapply(fmnl_settings$V, "-", maxV)
          fmnl_settings$V <- lapply(X=fmnl_settings$V, FUN=exp)
          fmnl_settings$V <- mapply('*', fmnl_settings$V, fmnl_settings$avail, SIMPLIFY = FALSE)
          denom = Reduce('+',fmnl_settings$V)
          P <- lapply(fmnl_settings$V, "/", denom)
        } else {
         ### work with subtracting the chosenV - USE THE ONE WITH THE HIGHEST SHARE
          chosenV <- mapply("*", fmnl_settings$maxShare, fmnl_settings$V, SIMPLIFY=FALSE)
          chosenV <- Reduce('+', chosenV)
          fmnl_settings$V <- lapply(X=fmnl_settings$V, "-", chosenV)
          fmnl_settings$V <- lapply(X=fmnl_settings$V, FUN=exp)
          # consider availabilities (it assumes V and avail are in the same order)
          fmnl_settings$V <- mapply('*', fmnl_settings$V, fmnl_settings$avail, SIMPLIFY = FALSE)
          # calculate the denominator of the Logit probability expression
          denom = Reduce('+',fmnl_settings$V)
          P <- lapply(fmnl_settings$V, "/", denom)
          if(any(sapply(P, anyNA))){
            P <- lapply(P, function(p){p[is.na(p)] <- 1; return(p)} )
            sP <- Reduce("+", P)
            P <- lapply(P, "/", sP)
          }
        }
        if(!fmnl_settings$choiceNA) P[["chosen"]] <- Reduce("*", mapply("^", P, fmnl_settings$choiceShares, SIMPLIFY=FALSE))
        # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
      } else { 
        ### work with subtracting the chosenV - USE THE ONE WITH THE HIGHEST SHARE
        chosenV <- mapply("*", fmnl_settings$maxShare, fmnl_settings$V, SIMPLIFY=FALSE)
        chosenV <- Reduce('+', chosenV)
        fmnl_settings$V <- lapply(X=fmnl_settings$V, "-", chosenV)
        numerator       <- exp(Reduce("+", mapply("*", fmnl_settings$V, fmnl_settings$choiceShares, SIMPLIFY=FALSE)))            
        fmnl_settings$V <- lapply(X=fmnl_settings$V, FUN=exp)
        # consider availabilities (it assumes V and avail are in the same order)
        fmnl_settings$V <- mapply('*', fmnl_settings$V, fmnl_settings$avail, SIMPLIFY = FALSE)
        # calculate the denominator of the Logit probability expression
        denom = Reduce('+',fmnl_settings$V)
        ### SH TO CHANGE
        P <- numerator/denom
      }
      ## insert excluded rows with value 1 if only teh chosen is requested, and 0 if all
      #if(any(!fmnl_settings$rows) & restoreRows){
      #  if(is.list(P)) P <- lapply(P, apollo_insertRows, r=fmnl_settings$r, val=0) else P <- apollo_insertRows(P, fmnl_settings$rows, 1)
      #}
      return(P)
    }
    
    fmnl_settings$fmnl_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
    
        # Table describing dependent variable
        choicematrix <- matrix(0, nrow=4, ncol=inputs$nAlt, 
                               dimnames=list(c("Times available",
                                               "Observations with non-zero share",
                                               "Average share overall",
                                               "Average share when available"),
                                             inputs$altnames))
        for(a in 1:inputs$nAlt){
          choicematrix[1,a] <- ifelse(length(inputs$avail[[a]])==1 && inputs$avail[[a]]==1, 
                                      inputs$nObs, sum(inputs$avail[[a]]) )
          choicematrix[2,a] <- sum(inputs$choiceShares[[a]]>0)
          choicematrix[3,a] <- ifelse(choicematrix[1,a]>0, mean(inputs$choiceShares[[a]]), 0)
          choicematrix[4,a] <- ifelse(choicematrix[1,a]>0, sum(inputs$choiceShares[[a]])/choicematrix[1,a], 0)
        }
        # Print table
        if(!apollo_inputs$silent & data){
          apollo_print("\n")
          apollo_print(paste0('Overview of choices for ', toupper(inputs$modeltype), ' model component ', 
                              ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
          print(round(choicematrix,2))
          
          # Print warnings
          for(a in 1:inputs$nAlt){
            if(choicematrix[2,a]==0) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'), type="w")
            if(choicematrix[4,a]==1) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always given the full choice share when available in model component "', inputs$componentName, '".'), type="w")
          }
          #if(inputs$avail_set==TRUE & !apollo_inputs$silent) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'), type="w")
        }
      
        return(invisible(TRUE))

    }
    
    
    # Store model type
    fmnl_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(fmnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    fmnl_settings$gradient <- FALSE
    if(test){
      fmnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, fmnl_settings$V)
      fmnl_settings$gradient <- !is.null(fmnl_settings$dV)
    }; rm(test)
    
    # Return fmnl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      fmnl_settings$V      <- NULL
      return(fmnl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(fmnl_settings$V, is.function))){
    fmnl_settings$V = lapply(fmnl_settings$V, function(f) if(is.function(f)) f() else f )
  } 
  fmnl_settings$V <- lapply(fmnl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V and drop rows if neccesary
  fmnl_settings$V <- fmnl_settings$V[fmnl_settings$altnames]
  if(!all(fmnl_settings$rows)) fmnl_settings$V <- lapply(fmnl_settings$V, apollo_keepRows, r=fmnl_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Check that alternatives are named in altcodes and V
    if(is.null(fmnl_settings$altnames) || is.null(fmnl_settings$altcodes) || is.null(names(fmnl_settings$V))) stop("SYNTAX ISSUE - Alternatives for model component \"",fmnl_settings$componentName,"\" must be named, both in 'alternatives' and 'utilities'.")
    
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(fmnl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) fmnl_settings$fmnl_diagnostics(fmnl_settings, apollo_inputs)
    
    testL = fmnl_settings$probs_fmnl(fmnl_settings, all=FALSE)
    if(any(!fmnl_settings$rows)) testL <- apollo_insertRows(testL, fmnl_settings$rows, 1)
    if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",fmnl_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",fmnl_settings$componentName,"\"", sep=""), type="i")
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    for(i in 1:length(fmnl_settings$avail)) if(length(fmnl_settings$avail[[i]])==1) fmnl_settings$avail[[i]] <- rep(fmnl_settings$avail[[i]], fmnl_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, fmnl_settings$avail)) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    if(any(!fmnl_settings$rows)) P <- apollo_insertRows(P, fmnl_settings$rows, 1)
    return(P)
  }
  
  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality=="shares_LL"){
    for(i in 1:length(mnl_settings$avail)) if(length(mnl_settings$avail[[i]])==1) mnl_settings$avail[[i]] <- rep(mnl_settings$avail[[i]], mnl_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, mnl_settings$avail)) # number of available alts in each observation
    Y = do.call(cbind,mnl_settings$Y)
    if(var(nAvAlt)==0){
      Yshares = colSums(Y)/nrow(Y)
      P = as.vector(Y%*%Yshares)
    } else {
      ## Estimate model with constants only
      mnl_ll = function(b, A, Y) as.vector(Y%*%c(b,0) - log(rowSums( A%*%exp(c(b,0)) )))
      A = do.call(cbind, mnl_settings$avail)
      b = maxLik::maxLik(mnl_ll, start=rep(0, mnl_settings$nAlt - 1), 
                         method='BFGS', finalHessian=FALSE, A=A, Y=Y)$estimate
      P = exp(mnl_ll(b, A, Y))
    }
    if(any(!mnl_settings$rows)) P <- apollo_insertRows(P, mnl_settings$rows, 1)
    return(P)
  }
  
  # ################################################################# #
  #### functionality="estimate/prediction/conditionals/raw/output" ####
  # ################################################################# #
  
  if(functionality %in% c("estimate","conditionals", "components", "output")){
    P <- fmnl_settings$probs_fmnl(fmnl_settings, all=FALSE)
    if(any(!fmnl_settings$rows)) P <- apollo_insertRows(P, fmnl_settings$rows, 1)
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- fmnl_settings$probs_fmnl(fmnl_settings, all=TRUE)
    if(any(!fmnl_settings$rows)) P <- lapply(P, apollo_insertRows, r=fmnl_settings$rows, val=NA)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
   # Verify everything necessary is available
   if(is.null(fmnl_settings$dV) || !all(sapply(fmnl_settings$dV, is.function))) stop("INTERNAL ISSUE - Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
   apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                           error=function(e) stop("INTERNAL ISSUE - apollo_fmnl could not fetch apollo_beta for gradient estimation."))
   if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_fmnl could not fetch apollo_inputs$database for gradient estimation.")
   
   # Calculate probabilities and derivatives of utilities for all alternatives
   P    <- fmnl_settings$probs_fmnl(fmnl_settings, all=TRUE)
   Pcho <- P[["chosen"]]
   P    <- P[-which(names(P)=="chosen")]
   e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
   for(i in 1:length(fmnl_settings$dV)) environment(fmnl_settings$dV[[i]]) <- e
   dV<- lapply(fmnl_settings$dV, function(dv) dv())
   if(!all(fmnl_settings$rows)) for(i in 1:length(dV)) dV[[i]] <- lapply(dV[[i]], apollo_keepRows, fmnl_settings$rows)
   for(i in 1:fmnl_settings$nAlt) dV[[i]] <- lapply(dV[[i]], 
                                                   function(dvik){ # Make dV=0 for unavailable alternatives
                                                     test <- length(dvik)==1 && length(fmnl_settings$avail[[i]])>1
                                                     if(test) dvik <- rep(dvik, fmnl_settings$nObs)
                                                     dvik[!fmnl_settings$avail[[i]]] <- 0
                                                     return(dvik)
                                                   })
   
   # Calculate gradient
   G_alts=list()
   for(s in 1:length(P)){
     GA <- list()
     for(j in 1:length(P)){
      if(j==s){
        GA[[j]]=1-P[[j]]
      } else {
        GA[[j]]=-P[[j]]
      }
     }
     G <- list()
     for(k in 1:length(dV[[1]])){
       dVk   <- lapply(dV, function(dv) dv[[k]])
       G[[k]] <- Reduce("+", mapply("*", GA, dVk, SIMPLIFY=FALSE))
     }
     G_alts[[s]] <- lapply(G, "*", fmnl_settings$choiceShares[[s]])
   }
   G=list()
   for(k in 1:length(dV[[1]])){
     G[[k]]=Pcho*Reduce("+",lapply(G_alts, function(g) g[[k]]))  
   }

   # Restore rows and return
   if(!all(fmnl_settings$rows)){
     Pcho <- apollo_insertRows(Pcho, fmnl_settings$rows, 1)
     G    <- lapply(G, apollo_insertRows, r=fmnl_settings$rows, val=0)
   }
   return(list(like=Pcho, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(fmnl_settings$fmnl_diagnostics(fmnl_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(fmnl_settings$fmnl_diagnostics(fmnl_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
