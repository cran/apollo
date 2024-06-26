#' Calculates the likelihood of a latent class model
#' 
#' Given within class probabilities, and class allocation probabilities, calculates the probabilities of an Exploded Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' 
#' @param lc_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                  \itemize{
#'                    \item \strong{inClassProb}: List of probabilities. Conditional likelihood for each class. One element per class, in the same order as \code{classProb}.
#'                    \item \strong{classProb}: List of probabilities. Allocation probability for each class. One element per class, in the same order as \code{inClassProb}.
#'                    \item \strong{componentName}: Character. Name given to model component (optional).
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
                          gradient      = TRUE,
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
  nRowsClassProb <- 0
  if(is.vector(classProb[[1]])) nRowsClassProb <- length(classProb[[1]]) else {
    if(functionality=="gradient"){
      nRowsClassProb <- dim(classProb[[1]][["like"]])[1]
      } else nRowsClassProb <- dim(classProb[[1]])[1]
  }
  
  # Count number of rows in inClassProb
  nRowsInClassProb <- 0
  ### SH bug fix 10 Oct - needed as a list is also a vector!
  #if(is.vector(inClassProb[[1]])) nRowsInClassProb <- length(inClassProb[[1]]) else {
  #  nRowsInClassProb <- dim(inClassProb[[1]])[1]
  #}
  ### DP bug fix 13 Aug 2021
  #if(is.list(inClassProb[[1]])){
  #  nRowsInClassProb <- length(inClassProb[[1]][[1]])
  #} else if(is.vector(inClassProb[[1]])) nRowsInClassProb <- length(inClassProb[[1]])
  if(is.list(inClassProb[[1]])) tmp <- inClassProb[[1]][[1]] else tmp <- inClassProb[[1]]
  if(is.array(tmp)) nRowsInClassProb <- dim(tmp)[1] else nRowsInClassProb <- length(tmp)
  rm(tmp)
  
  
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
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- setNames(sapply(as.list(unique(indivID)),function(x) sum(indivID==x)),unique(indivID))
        
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        classProb <- lapply(classProb, rowsum, group=indivID, reorder=FALSE)
        classProb <- lapply(classProb, function(x) if(is.matrix(x) && ncol(x)==1) return(as.vector(x)) else return(x))
        classProb <- lapply(classProb, '/', nObsPerIndiv)
        #if(!apollo_inputs$silent) apollo_print(paste0('Class probabilities for model component "', lc_settings$componentName, 
        #                                              '" averaged across observations of each individual.'))
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        classProb <- lapply(classProb, rep, nObsPerIndiv)
      }
    }
    
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
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb  & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- setNames(sapply(as.list(unique(indivID)),function(x) sum(indivID==x)),unique(indivID))
    
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        stop("SPECIFICATION ISSUE - Class-probability variable for model component \"",lc_settings$componentName,"\" has more elements than in-class-probability.")
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        S=length(classProb)
        for(s in 1:S){
          isMat <- is.matrix(classProb[[s]])
          isCub <- is.array(classProb[[s]]) && !isMat && length(dim(classProb[[s]]))==3
          if(isCub){
            # first average out over intra
            classProb[[s]]=colSums(aperm(classProb[[s]], perm=c(3,1,2)))/dim(classProb[[s]])[3]
          } 
          # now check what's left
          isVec <- is.vector(classProb[[s]])
          isMat <- is.matrix(classProb[[s]])
          if(isVec) classProb[[s]]=rep(classProb[[s]],nObsPerIndiv)
          if(isMat){
            tmp <- matrix(0, nrow=sum(nObsPerIndiv), ncol=ncol(classProb[[s]]))
            for(n in 1:length(nObsPerIndiv)){
              a <- ifelse(n==1, 1, sum(nObsPerIndiv[1:(n-1)]) + 1)
              b <- a + nObsPerIndiv[n] - 1 
              tmp[a:b,] <- rep(classProb[[s]][n,], each=nObsPerIndiv[n])
            }
            classProb[[s]] <- tmp
          } 
        }
        
      }
    }
    
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
    if(any(sapply(classProb, function(p) is.null(p$like) || is.null(p$grad)))) stop("INTERNAL ISSUE - Some components in classProb are missing the 'like' and/or 'grad' elements")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("INCORRECT FUNCTION/SETTING USE - workInLogs cannot be used in conjunction with analyticGrad")
    K <- length(inClassProb[[1]]$grad) # number of parameters
    if(K!=length(classProb[[1]]$grad)) stop("INTERNAL ISSUE - Dimensions of gradients not consistent between classProb and inClassProb")
    if(any(sapply(inClassProb, function(p) length(p$grad))!=K)) stop("INTERNAL ISSUE - Dimensions of gradients from different components in inClassProb imply different number of parameters")
    if(any(sapply(classProb, function(p) length(p$grad))!=K)) stop("INTERNAL ISSUE - Dimensions of gradients from different components in classProb imply different number of parameters")
    
    # Count number of rows in inClassProb
    nRowsInClassProb <- 0
    if(is.vector(inClassProb[[1]]$like)) nRowsInClassProb <- length(inClassProb[[1]]$like) else {
      nRowsInClassProb <- dim(inClassProb[[1]][["like"]])[1]
    }
    
    # Count number of obs in class alloc probabilities
    nRowsClassProb <- max(sapply(classProb, function(p) length(p$like)))
    # Match the number of rows in each list
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- setNames(sapply(as.list(unique(indivID)),function(x) sum(indivID==x)),unique(indivID))

      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        for(i in 1:length(classProb)){
          classProb[[i]]$grad <- lapply(classProb[[i]]$grad, rowsum, group=indivID, reorder=FALSE)
          classProb[[i]]$grad <- lapply(classProb[[i]]$grad, '/', nObsPerIndiv)
          classProb[[i]]$like <- rowsum(classProb[[i]]$like, group=indivID, reorder=FALSE)/nObsPerIndiv
        }
        #if(!apollo_inputs$silent) apollo_print(paste0('Class probabilities for model component "', lc_settings$componentName, 
        #                                              '" averaged across observations of each individual.'))
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        for(i in 1:length(classProb)){
          classProb[[i]]$grad <- lapply(classProb[[i]]$grad, rep, nObsPerIndiv)
          classProb[[i]]$like <- rep(classProb[[i]]$like, nObsPerIndiv)
        }
      }
    }
    
    ### Calculate actual likelihood and gradient
    nClass <- length(classProb)
    nParam <- length(inClassProb[[1]]$grad)
    Pout = list(like=0, grad=list())
    for(k in 1:nParam) Pout$grad[[k]] <- 0
    for(i in 1:nClass){
      Pout$like <- Pout$like + inClassProb[[i]]$like*classProb[[i]]$like
      for(k in 1:nParam) Pout$grad[[k]] <- Pout$grad[[k]] + 
          inClassProb[[i]]$grad[[k]]*classProb[[i]]$like + inClassProb[[i]]$like*classProb[[i]]$grad[[k]]
    }
    
    return(Pout)
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