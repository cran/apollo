#' Calculates class allocation probabilities for a Latent Class model
#'
#' Calculates class allocation probabilities for a Latent Class model using a Multinomial Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' @param classAlloc_settings List of inputs of the MNL model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the classes in class allocation model. Names of elements must match those in \code{avail}, if provided.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of classes, one element per class Names of elements must match those in \code{classes}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}. 
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                     }
#' @return The returned object depends on the value of argument \code{functionality}, which it fetches from the calling stack (see \link{apollo_validateInputs}).
#'         \itemize{
#'           \item \strong{\code{"components"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"estimate"}}: List of vector/matrices/arrays with the allocation probabilities for each class.
#'           \item \strong{\code{"gradient"}}: List containing the likelihood and gradient of the model component.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"prediction"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{classAlloc_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"report"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"shares_LL"}}: List with probabilities for each class in an equal shares setting.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: List with probabilities for each class in an equal shares setting.
#'         }
#' @export
#' @importFrom utils capture.output
apollo_classAlloc <- function(classAlloc_settings){
  # Fetch functionality
  functionality = tryCatch(get('functionality', parent.frame(n=2), inherits=TRUE), error=function(e) return('estimate'))
  
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE, silent=FALSE, analyticGrad=TRUE),
                                                          silent=FALSE) ))
  
  ### Set or extract componentName
  modelType = "classAlloc"
  if(is.null(classAlloc_settings[["componentName"]])){
    classAlloc_settings[["componentName"]] = ifelse(!is.null(classAlloc_settings[['componentName2']]),
                                             classAlloc_settings[['componentName2']], modelType)
    test <- functionality=="validate" && classAlloc_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', classAlloc_settings[["componentName"]], '" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, classAlloc_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", classAlloc_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used, and "alternatives" by "classes".
  if(!is.null(classAlloc_settings[["utilities"]])) names(classAlloc_settings)[which(names(classAlloc_settings)=="utilities")]="V"
  if(!is.null(classAlloc_settings[["alternatives"]])) names(classAlloc_settings)[which(names(classAlloc_settings)=="alternatives")]="classes"
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Search for any LC _settings that contain a classAlloc_settings element in it
  tmp <- grep('_settings$', ls(apollo_inputs), value=TRUE); tmp2 <- NULL
  for(i in tmp) if(!is.null(apollo_inputs[[i]][['classAlloc_settings']])) tmp2 <- c(tmp2, i)
  # Load the preprocessed settings only if there is a single one for calssAlloc
  if( length(tmp2)==1 && (functionality!="preprocess") ){
    # Load classAlloc_settings from apollo_inputs
    tmp <- apollo_inputs[[tmp2]][['classAlloc_settings']]
    # If there is no V inside the loaded classAlloc_settings, restore the one received as argument
    if(is.null(tmp$V)) tmp$V <- classAlloc_settings$V
    classAlloc_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    classAlloc_settings <- apollo_preprocess(inputs = classAlloc_settings, modelType, functionality, apollo_inputs)
    
    # Determine which mnl likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation for classAlloc available")
    # Using R likelihood
    classAlloc_settings$probs_MNL <- function(classAlloc_settings){
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      classAlloc_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                               classAlloc_settings$V, classAlloc_settings$avail, SIMPLIFY=FALSE)
      # subtract the maxV
      maxV <- do.call(pmax, classAlloc_settings$V)
      classAlloc_settings$V <- lapply(classAlloc_settings$V, "-", maxV)
      # calculate probabilities of all alternatives
      classAlloc_settings$V <- lapply(X=classAlloc_settings$V, FUN=exp)
      classAlloc_settings$V <- mapply('*', classAlloc_settings$V, classAlloc_settings$avail, SIMPLIFY = FALSE)
      denom = Reduce('+',classAlloc_settings$V)
      P <- lapply(classAlloc_settings$V, "/", denom)
      return(P)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(classAlloc_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    classAlloc_settings$gradient <- FALSE
    if(test){
      classAlloc_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, classAlloc_settings$V)
      classAlloc_settings$gradient <- !is.null(classAlloc_settings$dV)
    }; rm(test)
    
    # Return classAlloc_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      classAlloc_settings$V      <- NULL
      return(classAlloc_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(classAlloc_settings$V, is.function))){
    classAlloc_settings$V = lapply(classAlloc_settings$V, function(f) if(is.function(f)) f() else f )
  } 
  classAlloc_settings$V <- lapply(classAlloc_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Drop rows from V if necessary
  if(!all(classAlloc_settings$rows)) classAlloc_settings$V <- lapply(classAlloc_settings$V, apollo_keepRows, r=classAlloc_settings$rows)
  # No need to drop rows in avail, as it was already filtered durin pre-processing
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(classAlloc_settings, modelType, 
                                                                  functionality, apollo_inputs)
    
    # No diagnose for the class allocation component
    #if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(classAlloc_settings, modelType, apollo_inputs)
    
    testL = classAlloc_settings$probs_MNL(classAlloc_settings)
    #if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",classAlloc_settings$componentName,"\"")
    #if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",classAlloc_settings$componentName,"\"", sep=""))
    return(invisible(testL))
  }
  
  # ############################################ #
  #### functionality="zero_LL" or "shares_LL" ####
  # ############################################ #
  
  if(functionality %in% c("zero_LL", 'shares_LL')){
    # turn scalar availabilities into vectors
    for(i in 1:length(classAlloc_settings$avail)) if(length(classAlloc_settings$avail[[i]])==1) classAlloc_settings$avail[[i]] <- rep(classAlloc_settings$avail[[i]], classAlloc_settings$nObs) 
    P <- lapply(classAlloc_settings$avail, '/', Reduce('+', classAlloc_settings$avail))
    if(any(!classAlloc_settings$rows)) P <- lapply(P, apollo_insertRows, r=classAlloc_settings$rows, val=1)
    return(P)
  }
  
  # ################################################################################### #
  #### functionality="estimate/conditionals/components/output/prediction/raw/report" ####
  # ################################################################################### #
  
  if(functionality %in% c("estimate","conditionals", "components", "output", "prediction", "raw", "report")){
    P <- classAlloc_settings$probs_MNL(classAlloc_settings)
    if(!all(classAlloc_settings$rows)) P <- lapply(P, apollo_insertRows, r=classAlloc_settings$rows, val=1)
    #P <- lapply(P, apollo_firstRow, apollo_inputs=apollo_inputs)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necesary is available
    if(is.null(classAlloc_settings$dV) || !all(sapply(classAlloc_settings$dV, is.function))) stop("INTERNAL ISSUE - Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_mnl could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_mnl could not fetch apollo_inputs$database for gradient estimation.")
    # Calculate probabilities and derivatives of utilities for all alternatives
    P <- classAlloc_settings$probs_MNL(classAlloc_settings)
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    for(i in 1:length(classAlloc_settings$dV)) environment(classAlloc_settings$dV[[i]]) <- e
    dV<- lapply(classAlloc_settings$dV, function(dv) dv())
    if(!all(classAlloc_settings$rows)) for(i in 1:length(dV)) dV[[i]] <- lapply(dV[[i]], apollo_keepRows, classAlloc_settings$rows)
    for(i in 1:classAlloc_settings$nAlt) dV[[i]] <- lapply(dV[[i]], 
                                                    function(dvik){ # Make dV=0 for unavailable alternatives
                                                      test <- length(dvik)==1 && length(classAlloc_settings$avail[[i]])>1
                                                      if(test) dvik <- rep(dvik, classAlloc_settings$nObs)
                                                      dvik[!classAlloc_settings$avail[[i]]] <- 0
                                                      return(dvik)
                                                    })
    # Calculate gradient for each alternative
    K   <- length(dV[[1]])
    Pd <- list()
    for(k in 1:K) Pd[[k]] <- Reduce('+', mapply('*', P, lapply(dV, function(dv) dv[[k]]), SIMPLIFY=FALSE))
    G <- list()
    for(i in 1:classAlloc_settings$nAlt){
      G[[i]] <- mapply(function(dvik, pdk) dvik - pdk, dV[[i]], Pd, SIMPLIFY=FALSE)
      G[[i]] <- lapply(G[[i]],'*',P[[i]])
    }
    # Restore rows and return
    if(!all(classAlloc_settings$rows)){
      P <- lapply(P, apollo_insertRows, r=classAlloc_settings$rows, val=1)
      for(i in 1:classAlloc_settings$nAlt) G[[i]] <- lapply(G[[i]], apollo_insertRows, r=classAlloc_settings$rows, val=0)
    }
    # Change order and return
    ans <- list()
    for(i in 1:classAlloc_settings$nAlt) ans[[i]] <- list(like=P[[i]], grad=G[[i]])
    return(ans)
  }
  
}
