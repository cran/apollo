#' Calculates the likelihood of a latent class model
#' 
#' Given within class probabilities, and class allocation probabilities, calculates the probabilities of an Exploded Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' 
#' @param lc_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                  \itemize{
#'                    \item \strong{classProb}: List of probabilities. Allocation probability for each class. One element per class, in the same order as \code{inClassProb}.
#'                    \item \strong{componentName}: Character. Name given to model component (optional).
#'                    \item \strong{inClassProb}: List of probabilities. Conditional likelihood for each class. One element per class, in the same order as \code{classProb}.
#'                  }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
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
#'           \item \strong{\code{"components"}}: Returns nothing.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"gradient"}}: List containing the likelihood and gradient of the model component.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all models components, for each class.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{lc_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Class allocation overview.
#'           \item \strong{\code{"shares_LL"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but also runs a set of tests on the given arguments.
#'           \item \strong{\code{"zero_LL"}}: Same as \code{"estimate"}
#'         }
#' @importFrom utils capture.output
#' @export
apollo_lc <- function(lc_settings, apollo_inputs, functionality){
  modelType = 'LC'
  if(is.null(lc_settings[["componentName"]])){
    lc_settings[["componentName"]] = ifelse(!is.null(lc_settings[['componentName2']]),
                                            lc_settings[['componentName2']], modelType)
    test <- functionality=="validate" && lc_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', lc_settings[["componentName"]], '" by default.'))
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  if(!is.null(apollo_inputs[[paste0(lc_settings$componentName, "_settings")]]) && (functionality!="preprocess")){
    # Load lc_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(lc_settings$componentName, "_settings")]]
    # If there is no V inside the loaded mnl_settings, restore the one received as argument
    if(is.null(tmp$inClassProb)) tmp$inClassProb <- lc_settings$inClassProb
    if(is.null(tmp$classProb  )) tmp$classProb   <- lc_settings$classProb
    lc_settings <- tmp
    rm(tmp)
  } else {
    if(functionality=="preprocess"){
      # Do preprocessing
      countOutcomes <- function(L){ # Count number of outcomes inside EACH class
        if(!is.list(L)) return(0)
        if("rows" %in% names(L)) return(sum(L$rows))
        for(i in names(L)) return( sum(sapply(L, countOutcomes)) )
      }
      lc_settings <- list(componentName = lc_settings$componentName,
                          LCNObs        = max(sapply(lc_settings$inClassProb, countOutcomes)), 
                          LCCompNames   = names(lc_settings$inClassProb),
                          modelType     = 'LC',
                          gradient      = TRUE, # ADD CHECK
                          hessian       = TRUE, # ADD CHECK
                          classAlloc_settings = lc_settings$classProb)
    }
    # Diagnostics function
    lc_settings$lc_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      class_summary = matrix(0, nrow=length(inputs$classProb), ncol=1, dimnames=list(names(inputs$inClassProb), "Mean prob."))
      if(is.null(rownames(class_summary))) rownames(class_summary) <- paste0("Class_", 1:length(inputs$classProb))
      for(cc in 1:length(inputs$classProb)) class_summary[cc,1] <- mean(inputs$classProb[[cc]])
      if(!apollo_inputs$silent & param){
        apollo_print('\n')
        apollo_print(paste0('Summary of class allocation for ', toupper(inputs$modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        apollo_print(class_summary)
      }
      return(invisible(TRUE))
    }
    if(functionality=="preprocess") return(lc_settings)
  }
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality=="components") return(NULL)
  
  ### Special case for EM estimation
  test <- !is.null(apollo_inputs$EM) && apollo_inputs$EM 
  test <- test && !is.null(apollo_inputs$class_specific) && apollo_inputs$class_specific>0
  if(test) return(lc_settings$inClassProb[[1]])
  rm(test)

  # #################### #
  #### General checks ####
  # #################### #
  if(is.null(lc_settings[["componentName"]])){
    modelType <- 'LC'
    lc_settings[["componentName"]] = ifelse(!is.null(lc_settings[['componentName2']]),
                                            lc_settings[['componentName2']], modelType)
    test <- functionality=="validate" && lc_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 lc_settings[["componentName"]],'" by default.'))
  }
  if(is.null(lc_settings[["inClassProb"]])) stop('SYNTAX ISSUE - The lc_settings list for model component "',
                                                 lc_settings$componentName,'" needs to include an object called "inClassProb"!')
  if(is.null(lc_settings[["classProb"]])) stop('SYNTAX ISSUE - The lc_settings list for model component "',
                                               lc_settings$componentName,'" needs to include an object called "classProb"!')
  inClassProb    = lc_settings[["inClassProb"]]
  classProb      = lc_settings[["classProb"]]
  apollo_control = apollo_inputs[["apollo_control"]]
  if(apollo_control$workInLogs && apollo_control$mixing) stop('INTERNAL ISSUE - The settings "workInLogs" and "mixing" in "apollo_control" ',
                                                              'cannot be used together for latent class models.')
  
  # Check that inClassProb is not a list of lists (if its is, use 'model' component or fail)
  if(!is.list(inClassProb)) stop('INPUT ISSUE - Setting "inClassProb" inside "lc_settings" must be a list.')
  for(i in 1:length(inClassProb)) if(is.list(inClassProb[[i]]) && functionality %in% c("conditionals","estimate","validate","zero_LL", "shares_LL", "output")){
    test <- is.null(inClassProb[[i]]$model)
    if(test) stop(paste0('INPUT ISSUE - At least one element inside "inClassProb" setting is a list. It should be a numeric vector, matrix',
                         ' or array. If you are using apollo_combineModels inside each class, try using it with the argument', 
                         ' apollo_combineModel(..., asList=FALSE)')) else inClassProb[[i]] <- inClassProb[[i]]$model
  }
  
  # Count number of rows in classProb
  Dim1 <- function(x){if(is.array(x)) return(dim(x)[[1]]) else return(length(x))}
  nRowsClassProb <- Dim1(classProb[[1]])
  if(functionality %in% c("gradient", "hessian")) nRowsClassProb <- Dim1(classProb$like[[1]])
  
  # Count number of rows in inClassProb
  nRowsInClassProb <- Dim1(inClassProb[[1]])
  if(functionality %in% c("gradient", "hessian")) nRowsInClassProb <- Dim1(inClassProb[[1]]$like)
  
  # Count number of observations per individual
  indivID <- get(apollo_control$indivID)
  nObsPerIndiv <- unique(indivID)
  nObsPerIndiv <- setNames(sapply(as.list(nObsPerIndiv),
                                  function(x) sum(indivID==x)), nObsPerIndiv)
  
  # Function to resize classProbs
  resize <- function(x, newL){
    # If it is a list, apply function recursively
    if(is.list(x)) return( lapply(x, resize, newL=newL) )
    # If no change is needed, return as is
    if(is.array(x) && length(x)==1) x <- as.vector(x)
    if(is.array(x)) oldL <- dim(x)[1] else oldL <- length(x)
    isCube <- is.array(x) && length(dim(x))==3
    if(is.matrix(x) && ncol(x)==1) x <- as.vector(x)
    if(oldL==1 | oldL==newL) return(x)
    # If longer
    if(oldL>newL){
      if(isCube){
        x <- apply(x, MARGIN=c(2,3), FUN=rowsum, group=c(1,1,2), reorder=FALSE)
      } else x <- rowsum(x, group=indivID, reorder=FALSE)
      x <- x/nObsPerIndiv
    } 
    # If shorter
    if(oldL<newL){
      if(isCube){
        tmp <- array(0, dim=c(newL, dim(x)[2:3]))
        tmp[1:nObsPerIndiv[1],,] <- x[1,,]
        if(length(nObsPerIndiv)>1) for(n in 2:length(nObsPerIndiv)){
          tmp[sum(nObsPerIndiv[1:(n-1)]):sum(nObsPerIndiv[1:n]),,] <- x[n,,]
        }; x <- tmp; rm(tmp)
      } else if(is.matrix(x)){
        tmp <- matrix(0, newL, ncol(x))
        tmp[1:nObsPerIndiv[1],] <- x[1,]
        if(length(nObsPerIndiv)>1) for(n in 2:length(nObsPerIndiv)){
          tmp[sum(nObsPerIndiv[1:(n-1)]):sum(nObsPerIndiv[1:n]),] <- x[n,]
        }; x <- tmp; rm(tmp)
      } else x <- rep(x, times=nObsPerIndiv)
    }
    # Simplify x to a vector if appropriate
    if(is.matrix(x) && ncol(x)==1) x <- as.vector(x)
    return(x)
  }
  
  # ############## #
  #### Validate ####
  # ############## #
  if(functionality=="validate"){
    
    # Validate inputs
    if(!apollo_inputs$apollo_control$noValidation){
      if(length(inClassProb)!=length(classProb)) stop("INPUT ISSUE - Arguments 'inClassProb' and 'classProb' for model component \"",
                                                      lc_settings$componentName,"\" must have the same length.")
      for(cc in 1:length(classProb)) if(sum(inClassProb[[cc]])==0) stop('CALCULATION ISSUE - Within class probabilities for model component "',
                                                                        lc_settings$componentName, 
                                                                        '" are zero for every individual for class ',
                                                                        ifelse(is.numeric(names(inClassProb)[cc]),cc,names(inClassProb)[cc]))
      classprobsum = Reduce("+",classProb)
      if(any(round(classprobsum,10)!=1)) stop("SPECIFICATION ISSUE - Class allocation probabilities in 'classProb' for model component \"",
                                              lc_settings$componentName,"\" must sum to 1.")
      # check that probs are different across classes
      if(length(unique(sapply(inClassProb,sum)))!=length(inClassProb))  stop("INPUT ISSUE - At your starting values, the probabilities are the same across some of the classes.",
                                                                             "Please use different starting values across classes. If you still wish to run your model in the current form, please set noValidation=TRUE in apollo_control.")
    }
    
    
    # Print diagnostics
    if(!apollo_inputs$apollo_control$noDiagnostics) lc_settings$lc_diagnostics(lc_settings, apollo_inputs)
    
    # Check if classProb needs to be averaged to match the dimensionality of inClassProb, warn if so
    test <- nRowsInClassProb!=nRowsClassProb && nRowsClassProb!=1 && nRowsInClassProb < nRowsClassProb && !apollo_inputs$silent
    #if(test) apollo_print(paste0('Class probabilities for model component "', lc_settings$componentName, 
    #                             '" averaged across observations of each individual.'))
    if(test) apollo_print(paste0('The class allocation probabilities for model component "', lc_settings$componentName, 
                                 '" are calculated at the observation level in \'apollo_lcPars\', but are used in',
                                 ' \'apollo_probabilities\' to multiply within class probabilities that are at the',
                                 ' individual level. Apollo will average the class allocation probabilities across',
                                 ' observations for the same individual level before using them to multiply the',
                                 ' within-class probabilities. If your class allocation probabilities are',
                                 ' constant across choice situations for the same individual, then this is of no concern.',
                                 ' If your class allocation probabilities however vary across choice tasks, then you should',
                                 ' change your model specification in \'apollo_probabilities\' to only call \'apollo_panelProd\'',
                                 ' after calling \'apollo_lc\'.'))
    
    testL = apollo_lc(lc_settings, apollo_inputs, functionality="estimate")
    return(testL)
  } 
  
  
  # ######################################################## #
  #### Estimate, conditionals, output, zero_LL, shares_LL ####
  # ######################################################## #
  if(functionality %in% c("estimate", "conditionals", "output", "zero_LL", "shares_LL")){
    # Match the number of rows in each list
    classProb <- resize(classProb, nRowsInClassProb)
    # # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    # if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
    #     
    #   if( nRowsInClassProb < nRowsClassProb ){
    #     # If classProb is calculated at the obs level while the inClassProb is at the indiv level
    #     classProb <- lapply(classProb, rowsum, group=indivID, reorder=FALSE)
    #     classProb <- lapply(classProb, function(x) if(is.matrix(x) && ncol(x)==1) return(as.vector(x)) else return(x))
    #     classProb <- lapply(classProb, '/', nObsPerIndiv)
    #     #if(!apollo_inputs$silent) apollo_print(paste0('Class probabilities for model component "', lc_settings$componentName, 
    #     #                                              '" averaged across observations of each individual.'))
    #   } else {
    #     # If inClassProb is calculated at the obs level while the classProb is at the indiv level
    #     classProb <- lapply(classProb, rep, nObsPerIndiv)
    #   }
    # }
    
    # Calculate inClassProb*classProb
    nIndiv <- length(unique(get(apollo_control$indivID)))                                      ## FIX DP 18/05/2020
    checkForPanel <- function(x) (is.vector(x) && length(x)==nIndiv) || (is.array(x) && dim(x)[1]==nIndiv) ## FIX DP 13/08/2021
    panelProdApplied <- sapply(inClassProb, checkForPanel)                                                 ## FIX DP 13/08/2021
    if(apollo_control$workInLogs && all(panelProdApplied)){                                    ## FIX DP 18/05/2020
      # Assuming panelProd was applied inside each class, inClassProb is log(P[ns])            ## FIX DP 18/05/2020
      Bbar <- Reduce("+", inClassProb)/length(inClassProb)
      Pout <- mapply(function(pi, lnP) pi*exp(lnP - Bbar), classProb, inClassProb, SIMPLIFY=FALSE)
      Pout <- Bbar + log( Reduce("+", Pout) )
    } else {
      Pout <- mapply(function(p,pi) pi*p, inClassProb, classProb, SIMPLIFY=FALSE)
      Pout <- Reduce('+', Pout)
    }
    return( Pout )
  }
  
  # ##################### #
  #### Prediction, raw ####
  # ##################### #
  if(functionality %in% c("prediction", "raw")){
    # Match the number of rows in each list
    classProb <- resize(classProb, nRowsInClassProb)
    # # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    # if( nRowsInClassProb!=nRowsClassProb  & nRowsClassProb!=1 ){
    #   
    #   if( nRowsInClassProb < nRowsClassProb ){
    #     # If classProb is calculated at the obs level while the inClassProb is at the indiv level
    #     stop("SPECIFICATION ISSUE - Class-probability variable for model component \"",lc_settings$componentName,"\" has more elements than in-class-probability.")
    #   } else {
    #     # If inClassProb is calculated at the obs level while the classProb is at the indiv level
    #     S=length(classProb)
    #     for(s in 1:S){
    #       isMat <- is.matrix(classProb[[s]])
    #       isCub <- is.array(classProb[[s]]) && !isMat && length(dim(classProb[[s]]))==3
    #       if(isCub){
    #         # first average out over intra
    #         classProb[[s]]=colSums(aperm(classProb[[s]], perm=c(3,1,2)))/dim(classProb[[s]])[3]
    #       } 
    #       # now check what's left
    #       isVec <- is.vector(classProb[[s]])
    #       isMat <- is.matrix(classProb[[s]])
    #       if(isVec) classProb[[s]]=rep(classProb[[s]],nObsPerIndiv)
    #       if(isMat){
    #         tmp <- matrix(0, nrow=sum(nObsPerIndiv), ncol=ncol(classProb[[s]]))
    #         for(n in 1:length(nObsPerIndiv)){
    #           a <- ifelse(n==1, 1, sum(nObsPerIndiv[1:(n-1)]) + 1)
    #           b <- a + nObsPerIndiv[n] - 1 
    #           tmp[a:b,] <- rep(classProb[[s]][n,], each=nObsPerIndiv[n])
    #         }
    #         classProb[[s]] <- tmp
    #       } 
    #     }
    #     
    #   }
    # }
    
    nClass<- length(classProb)
    for(s in 1:length(inClassProb)){
      if(!is.list(inClassProb[[s]])) inClassProb[[s]]=list(inClassProb[[s]])  
    }
    
    nAlts <- length(inClassProb[[1]])
    
    Pout <- vector(mode="list", length=nAlts)
    for(i in 1:nAlts){
      Pout[[i]] <- 0*inClassProb[[1]][[1]]
      for(k in 1:nClass) Pout[[i]] <- Pout[[i]] + inClassProb[[k]][[i]]*classProb[[k]]
    }
    names(Pout)=names(inClassProb[[1]])
    return(Pout)
    
  }
  
  
  # ############## #
  #### Gradient ####
  # ############## #
  
  if(functionality==("gradient")){
    # Initial checks
    if(any(sapply(inClassProb, function(p) is.null(p$like) || is.null(p$grad)))) stop("INTERNAL ISSUE - Some components in inClassProb are missing the 'like' and/or 'grad' elements")
    if(is.null(classProb$like) || is.null(classProb$grad)) stop("INTERNAL ISSUE - Some components in classProb are missing the 'like' and/or 'grad' elements")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("INCORRECT FUNCTION/SETTING USE - workInLogs cannot be used in conjunction with analyticGrad")
    K <- length(inClassProb[[1]]$grad) # number of parameters
    if(length(classProb$grad)!=K) stop("INTERNAL ISSUE - Dimensions of gradients not consistent between classProb and inClassProb")
    if(any(sapply(inClassProb, function(p) length(p$grad))!=K)) stop("INTERNAL ISSUE - Dimensions of gradients from different components in inClassProb imply different number of parameters")
    if(length(classProb$grad)!=K) stop("INTERNAL ISSUE - Dimensions of gradients from different components in classProb imply different number of parameters")
    S <- length(classProb$like) # number of classes
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - Gradient ",
                                                   "could not fetch apollo_beta"))
    parName <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_inputs$apollo_fixed)]
    
    # Match the number of rows in each list
    classProb <- resize(classProb, nRowsInClassProb)
    inClassProb <- resize(inClassProb, nRowsInClassProb)
    # # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    # if( nRowsClassProb!=1 && nRowsInClassProb < nRowsClassProb){
    #   classProb$like <- lapply(classProb$like, \(l) rowsum(l, group=indivID, reorder=FALSE)/nObsPerIndiv)
    #   for(k in 1:K) classProb$grad[[k]] <- lapply(classProb$grad[[k]], \(g) rowsum(g, group=indivID, reorder=FALSE)/nObsPerIndiv)
    # }
    # if( nRowsClassProb!=1 && nRowsInClassProb > nRowsClassProb){
    #   classProb$like <- lapply(classProb$like, rep, times=nObsPerIndiv)
    #   for(k in 1:K) classProb$grad[[k]] <- lapply(classProb$grad[[k]], rep, nObsPerIndiv)
    # }
    
    ### Calculate actual likelihood and gradient
    nClass <- length(inClassProb)
    Pout = list(like=0, grad=setNames(vector(mode="list", length=K), parName))
    for(k in 1:K) Pout$grad[[k]] <- 0
    for(i in 1:nClass){
      Pout$like <- Pout$like + inClassProb[[i]]$like*classProb$like[[i]]
      for(k in 1:K) Pout$grad[[k]] <- Pout$grad[[k]] + 
          inClassProb[[i]]$grad[[k]]*classProb$like[[i]] + inClassProb[[i]]$like*classProb$grad[[k]][[i]]
      if(is.array(Pout$grad[[k]])) rownames(Pout$grad[[k]]) <- NULL
      if(is.vector(Pout$grad[[k]])) names(Pout$grad[[k]]) <- NULL
    }
    
    return(Pout)
  }
  
  # ############# #
  #### Hessian ####
  # ############# #
  
  if(functionality==("hessian")){
    ### TO DO
     # Checks
     # No need to check for rows
    K <- length(inClassProb[[1]]$grad) # number of parameters
    S <- length(classProb$like) # number of classes
    
    # Match the number of rows in each list
    classProb <- resize(classProb, nRowsInClassProb)
    # # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    # if( nRowsClassProb!=1 && nRowsInClassProb < nRowsClassProb){
    #   classProb$like <- lapply(classProb$like, \(l) rowsum(l, group=indivID, reorder=FALSE)/nObsPerIndiv)
    #   for(k in 1:K){
    #     classProb$grad[[k]] <- lapply(classProb$grad[[k]], \(g) rowsum(g, group=indivID, reorder=FALSE)/nObsPerIndiv)
    #     for(k2 in 1:K) classProb$hess[[k]][[k2]] <- lapply(classProb$hess[[k]][[k2]], \(h) rowsum(h, group=indivID, reorder=FALSE)/nObsPerIndiv)
    #   }; rm(k, k2)
    # }
    # if( nRowsClassProb!=1 && nRowsInClassProb > nRowsClassProb){
    #   classProb$like <- lapply(classProb$like, rep, times=nObsPerIndiv)
    #   for(k in 1:K){
    #     classProb$grad[[k]] <- lapply(classProb$grad[[k]], rep, times=nObsPerIndiv)
    #     for(k2 in 1:K) classProb$hess[[k]][[k2]] <- lapply(classProb$hess[[k]][[k2]], rep, times=nObsPerIndiv)
    #   }; rm(k, k2)
    # }
    
    # Calculate likelihood and gradient
    nClass <- length(classProb$like)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - Gradient ",
                                                   "could not fetch apollo_beta"))
    parName <- names(apollo_beta)[!(names(apollo_beta) %in% apollo_inputs$apollo_fixed)]
    L <- 0
    G <- setNames(vector(mode="list", length=K), parName)
    for(k in 1:K) G[[k]] <- 0
    for(s in 1:S){
      L <- L + inClassProb[[s]]$like*classProb$like[[s]]
      if(is.array(L)) rownames(L) <- NULL else names(L) <- NULL
      for(k in 1:K) G[[k]] <- G[[k]] + 
          inClassProb[[s]]$grad[[k]]*classProb$like[[s]] + inClassProb[[s]]$like*classProb$grad[[k]][[s]]
      if(is.array(G[[k]])) rownames(G[[k]]) <- NULL else names(G[[k]]) <- NULL
    }
    
    # Calculate the hessian
    H <- setNames(vector(mode="list", length=K), parName)
    for(k1 in 1:K){
      H[[k1]] <- setNames(vector(mode="list", length=K), parName)
      for(k2 in 1:K){
        H[[k1]][[k2]] <- 0
        for(s in 1:S) H[[k1]][[k2]] <- H[[k1]][[k2]] + 
            classProb$hess[[k1]][[k2]][[s]]*inClassProb[[s]]$like + 
            classProb$grad[[k1]][[s]]*inClassProb[[s]]$grad[[k2]] + 
            classProb$grad[[k2]][[s]]*inClassProb[[s]]$grad[[k1]] + 
            classProb$like[[s]]*inClassProb[[s]]$hess[[k1]][[k2]]
        if(is.array(H[[k1]][[k2]])) rownames(H[[k1]][[k2]]) <- NULL else 
          names(H[[k1]][[k2]]) <- NULL
      }
    }
    
    return(list(like=L, grad=G, hess=H))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(lc_settings$lc_diagnostics(lc_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(lc_settings$lc_diagnostics(lc_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}