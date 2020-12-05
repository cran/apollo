#' Calculates the likelihood of a latent class model
#' 
#' Using the conditional likelihoods of each latent class, as well as teir classification probabilities, calculate the weighted likelihood of the whole model.
#' 
#' @param lc_settings List of arguments used by \code{apollo_lc}. It must include the following.
#'                  \itemize{
#'                    \item \strong{inClassProb}: List of probabilities. Conditional likelihood for each class. One element per class, in the same order as \code{classProb}.
#'                    \item \strong{classProb}: List of probabilities. Allocation probability for each class. One element per class, in the same order as \code{inClassProb}.
#'                    \item \strong{componentName}: Character. Name given to model component.
#'                  }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate" Used for model estimation.
#'                        \item "prediction" Used for model predictions.
#'                        \item "validate" Used for validating input.
#'                        \item "zero_LL" Used for calculating null likelihood.
#'                        \item "conditionals" Used for calculating conditionals.
#'                        \item "output" Used for preparing output after model estimation.
#'                        \item "raw" Used for debugging.
#'                        \item "components" Returns \code{P} without changes.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all models components, for each class.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but also runs a set of tests on the given arguments.
#'           \item \strong{\code{"zero_LL"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
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
  
  # ############################################### #
  #### functionalities with untransformed return ####
  # ############################################### #
  
  if(functionality=="components") return(NULL)
  
  if(functionality=="preprocess"){
    lc_settings <- list(componentName = lc_settings$componentName,
                        LCNObs        = max(sapply(lc_settings$inClassProb, function(m_settings) sum(m_settings$rows))), 
                        LCCompNames   = names(lc_settings$inClassProb))
    return(lc_settings)
  }
  
  ### Special case for EM estimation
  test <- !is.null(apollo_inputs$EM) && apollo_inputs$EM 
  test <- test && !is.null(apollo_inputs$class_specific) && apollo_inputs$class_specific>0
  if(test) return(lc_settings$inClassProb[[1]])
  rm(test)

  # #################### #
  #### General checks ####
  # #################### #
  if(is.null(lc_settings[["componentName"]])){
    modelType <- 'LatentClass'
    lc_settings[["componentName"]] = ifelse(!is.null(lc_settings[['componentName2']]),
                                            lc_settings[['componentName2']], modelType)
    test <- functionality=="validate" && lc_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 lc_settings[["componentName"]],'" by default.'))
  }
  if(is.null(lc_settings[["inClassProb"]])) stop('The lc_settings list for model component "',
                                                 lc_settings$componentName,'" needs to include an object called "inClassProb"!')
  if(is.null(lc_settings[["classProb"]])) stop('The lc_settings list for model component "',
                                               lc_settings$componentName,'" needs to include an object called "classProb"!')
  inClassProb    = lc_settings[["inClassProb"]]
  classProb      = lc_settings[["classProb"]]
  apollo_control = apollo_inputs[["apollo_control"]]
  if(apollo_control$workInLogs && apollo_control$mixing) stop('The settings "workInLogs" and "mixing" in "apollo_control" ',
                                                              'cannot be used together for latent class models.')
  
  # Check that inClassProb is not a list of lists (if its is, use 'model' component or fail)
  if(!is.list(inClassProb)) stop('Setting "inClassProb" inside "lc_settings" must be a list.')
  for(i in 1:length(inClassProb)) if(is.list(inClassProb[[i]]) && functionality %in% c("conditionals","estimate","validate","zero_LL", "output")){
    test <- is.null(inClassProb[[i]]$model)
    if(test) stop(paste0('At least one element inside "inClassProb" setting is a list. It should be a numeric vector, matrix',
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
  if(is.list(inClassProb[[1]])){
    nRowsInClassProb <- length(inClassProb[[1]][[1]])
  } else if(is.vector(inClassProb[[1]])) nRowsInClassProb <- length(inClassProb[[1]])
  
  # ############## #
  #### Validate ####
  # ############## #
  if(functionality=="validate"){
    
    if(!apollo_inputs$apollo_control$noValidation){
      if(length(inClassProb)!=length(classProb)) stop("Arguments 'inClassProb' and 'classProb' for model component \"",
                                                      lc_settings$componentName,"\" must have the same length.")
      for(cc in 1:length(classProb)) if(sum(inClassProb[[cc]])==0) stop('Within class probabilities for model component "',
                                                                        lc_settings$componentName, 
                                                                        '" are zero for every individual for class ',
                                                                        ifelse(is.numeric(names(inClassProb)[cc]),cc,names(inClassProb)[cc]))
      classprobsum = Reduce("+",classProb)
      if(any(round(classprobsum,10)!=1)) stop("Class allocation probabilities in 'classProb' for model component \"",
                                              lc_settings$componentName,"\" must sum to 1.")
    }
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(lc_settings, modelType="LC", apollo_inputs)
    
    testL = apollo_lc(lc_settings, apollo_inputs, functionality="estimate")
    return(testL)
  } 
  
  
  # ############################################# #
  #### Estimate, conditionals, output, zero_LL ####
  # ############################################# #
  if(functionality %in% c("estimate", "conditionals", "output", "zero_LL")){
    # Match the number of rows in each list
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- as.vector(table(indivID))
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        classProb <- lapply(classProb, rowsum, group=indivID)
        classProb <- lapply(classProb, '/', nObsPerIndiv)
        if(!apollo_inputs$silent) apollo_print(paste0('Class probability for model component "', 
                                                      lc_settings$componentName, 
                                                      '" averaged across observations of each individual.'))
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        classProb <- lapply(classProb, rep, nObsPerIndiv)
        ### warning('Class probability repeated for all observations of the same individual.')
      }
    }
    
    # Calculate inClassProb*classProb
    nIndiv <- length(unique(get(apollo_control$indivID)))                                  ## FIX DP 18/05/2020
    panelProdApplied <- sapply(inClassProb, function(v) is.vector(v) && length(v)==nIndiv) ## FIX DP 18/05/2020
    if(apollo_control$workInLogs && all(panelProdApplied)){                                ## FIX DP 18/05/2020
      # Assuming panelProd was applied inside each class, inClasProb is log(P[ns])         ## FIX DP 18/05/2020
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
      nObsPerIndiv <- as.vector(table(indivID))
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        stop("Class-probability variable for model component \"",lc_settings$componentName,"\" has more elements than in-class-probability.")
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
    if(any(sapply(inClassProb, function(p) is.null(p$like) || is.null(p$grad)))) stop("Some components in inClassProb are missing the 'like' and/or 'grad' elements")
    if(any(sapply(classProb, function(p) is.null(p$like) || is.null(p$grad)))) stop("Some components in classProb are missing the 'like' and/or 'grad' elements")
    if(apollo_inputs$apollo_control$workInLogs && apollo_inputs$apollo_control$analyticGrad) stop("workInLogs cannot be used in conjunction with analyticGrad")
    K <- length(inClassProb[[1]]$grad) # number of parameters
    if(K!=length(classProb[[1]]$grad)) stop("Dimensions of gradients not consistent between classProb and inClassProb")
    if(any(sapply(inClassProb, function(p) length(inClassProb$grad))!=K)) stop("Dimensions of gradients from different components in inClassProb imply different number of parameters")
    if(any(sapply(classProb, function(p) length(classProb$grad))!=K)) stop("Dimensions of gradients from different components in classProb imply different number of parameters")
    
    # Count number of rows in inClassProb
    nRowsInClassProb <- 0
    if(is.vector(inClassProb[[1]]$like)) nRowsInClassProb <- length(inClassProb[[1]]$like) else {
      nRowsInClassProb <- dim(inClassProb[[1]][["like"]])[1]
    }
    
    # Match the number of rows in each list
    # The dimension of classProb is changed if necessary, dim(inClassProb) remains the same
    if( nRowsInClassProb!=nRowsClassProb & nRowsClassProb!=1 ){
      indivID <- get(apollo_control$indivID) 
      nObsPerIndiv <- as.vector(table(indivID))
      if( nRowsInClassProb < nRowsClassProb ){
        # If classProb is calculated at the obs level while the inClassProb is at the indiv level
        classProb$grad <- lapply(classProb$grad, rowsum, group=indivID)
        classProb$grad <- lapply(classProb$grad, '/', nObsPerIndiv)
        classProb$like <- lapply(classProb$like, rowsum, group=indivID)
        classProb$like <- lapply(classProb$like, '/', nObsPerIndiv)
        if(!apollo_inputs$silent) apollo_print(paste0('Class probability for model component "', 
                                                      lc_settings$componentName, 
                                                      '" averaged across observations of each individual.'))
      } else {
        # If inClassProb is calculated at the obs level while the classProb is at the indiv level
        classProb$grad <- lapply(classProb$grad, rep, nObsPerIndiv)
        classProb$like <- lapply(classProb$like, rep, nObsPerIndiv)
        ### warning('Class probability repeated for all observations of the same individual.')
      }
    }
    
    Pout = list()
    Pout$like=Reduce('+',mapply(function(p,pi) pi*p, inClassProb$like, classProb$like, SIMPLIFY=FALSE))
    Pout$like=Reduce('+',mapply(function(plike, pilike ,pgrad ,pigrad) pigrad*plike + pilike*pgrad, 
                                inClassProb$like, classProb$like, inClassProb$grad, classProb$grad, SIMPLIFY=FALSE))

    return(Pout)
    
  } 
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(lc_settings, modelType="LC", apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(lc_settings, modelType="LC", apollo_inputs, data =FALSE))
    return(P)
  }
  
}