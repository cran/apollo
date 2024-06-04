#' Calculates density for a Normal distribution
#'
#' Calculates density for a Normal distribution at a specific value with a specified mean and standard deviation and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' This function calcualtes the probability of the linear model outcomeNormal = mu + xNormal + epsilon, where epsilon is a random error distributed Normal(0,sigma).
#' If using this function in the context of an Integrated Choice and Latent Variable (ICLV) model with continuous
#' indicators, then \code{outcomeNormal} would be the value of the indicator, \code{xNormal} would be the value of the latent variable (possibly
#' multiplied by a parameter to measure its correlation with the indicator, e.g. xNormal=lambda*LV), and \code{mu} would be
#' an additional parameter to be estimated (the mean of the indicator, which should be fixed to zero if the indicator is
#' centered around its mean beforehand).
#' @param normalDensity_settings List of arguments to the functions. It must contain the following.
#'                               \itemize{
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                                 \item \strong{\code{mu}}: Numeric scalar. Intercept of the linear model.
#'                                 \item \strong{\code{outcomeNormal}}: Numeric vector. Dependent variable.
#'                                 \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                 \item \strong{\code{sigma}}: Numeric scalar. Variance of error component of linear model to be estimated.
#'                                 \item \strong{\code{xNormal}}: Numeric vector. Single explanatory variable.
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
#'           \item \strong{\code{"gradient"}}: Not implemented
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: Not implemented. Returns NA.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{normalDensity_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"report"}}: Dependent variable overview.
#'           \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'         }
#' @importFrom utils capture.output
#' @export
apollo_normalDensity <- function(normalDensity_settings, functionality){
  ### Set or extract componentName
  modelType = "NormD"
  if(is.null(normalDensity_settings[["componentName"]])){
    normalDensity_settings[["componentName"]] = ifelse(!is.null(normalDensity_settings[['componentName2']]),
                                                       normalDensity_settings[['componentName2']], modelType)
    test <- functionality=="validate" && normalDensity_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 normalDensity_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, normalDensity_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", normalDensity_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  if( !is.null(apollo_inputs[[paste0(normalDensity_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(normalDensity_settings$componentName, "_settings")]]
    # Restore variable lements if they do not exist
    if(is.null(tmp$xNormal)) tmp$xNormal <- normalDensity_settings$xNormal
    if(is.null(tmp$mu     )) tmp$mu      <- normalDensity_settings$mu     
    if(is.null(tmp$sigma  )) tmp$sigma   <- normalDensity_settings$sigma  
    normalDensity_settings <- tmp
    rm(tmp)
  } else {
    # Do pre-processing common to most models
    normalDensity_settings = apollo_preprocess(normalDensity_settings, modelType, 
                                               functionality, apollo_inputs)
    # diagnostics function
    normalDensity_settings$normalDensity_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      if(!apollo_inputs$silent & data){
        apollo_print('\n')
        apollo_print(paste0('Summary statistics for ', toupper(inputs$modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        tmp <- t(as.matrix(summary(inputs$outcomeNormal))); rownames(tmp) <- ""
        print(round(tmp,4))
        rm(tmp)
      }
      return(invisible(TRUE))
    }
    # Store model type
    normalDensity_settings$modelType <- modelType
    # Construct necessary input for gradient
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && apollo_inputs$apollo_control$analyticGrad 
    test <- test && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.function(normalDensity_settings$xNormal) 
    test <- test && is.function(normalDensity_settings$mu)
    test <- test && is.function(normalDensity_settings$sigma)
    normalDensity_settings$gradient <- FALSE
    if(test){
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(normalDensity_settings$xNormal))
      normalDensity_settings$dX     <- lapply(tmp, `[[`, 1)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(normalDensity_settings$mu     ))
      normalDensity_settings$dMu    <- lapply(tmp, `[[`, 1)
      tmp <- apollo_dVdB(apollo_beta, apollo_inputs, list(normalDensity_settings$sigma  ))
      normalDensity_settings$dSigma <- lapply(tmp, `[[`, 1)
      tmp <-        all( sapply(normalDensity_settings$dX    , is.function) )
      tmp <- tmp && all( sapply(normalDensity_settings$dMu   , is.function) )
      tmp <- tmp && all( sapply(normalDensity_settings$dSigma, is.function) )
      normalDensity_settings$gradient <- tmp
    }; rm(test)
    
    # Return normalDensity_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      normalDensity_settings$xNormal <- NULL
      normalDensity_settings$mu      <- NULL
      normalDensity_settings$sigma   <- NULL
      return(normalDensity_settings)
    }
  }
  
  
  # ################################################### #
  #### Evaluate xNormal, mu & sigma, and remove rows ####
  # ################################################### #
  
  ### Evaluate xNormal, mu & sigma
  if(is.function(normalDensity_settings$xNormal)) normalDensity_settings$xNormal <- normalDensity_settings$xNormal()
  if(is.function(normalDensity_settings$mu     )) normalDensity_settings$mu      <- normalDensity_settings$mu()
  if(is.function(normalDensity_settings$sigma  )) normalDensity_settings$sigma   <- normalDensity_settings$sigma()
  
  ### Drop excluded rows
  if(any(!normalDensity_settings$rows)){
    normalDensity_settings$outcomeNormal <- apollo_keepRows(normalDensity_settings$outcomeNormal, normalDensity_settings$rows)
    normalDensity_settings$xNormal       <- apollo_keepRows(normalDensity_settings$xNormal      , normalDensity_settings$rows)
    normalDensity_settings$mu            <- apollo_keepRows(normalDensity_settings$mu           , normalDensity_settings$rows)
    normalDensity_settings$sigma         <- apollo_keepRows(normalDensity_settings$sigma        , normalDensity_settings$rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){

    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
                                error=function(e) return(list(noValidation=FALSE, noDiagnostics=FALSE)) )

    if(!apollo_control$noValidation) apollo_validate(normalDensity_settings, modelType,
                                                     functionality, apollo_inputs)
    if(!apollo_control$noDiagnostics) normalDensity_settings$normalDensity_diagnostics(normalDensity_settings, apollo_inputs)
    testL = stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                         normalDensity_settings$mu, normalDensity_settings$sigma)
    if(any(!normalDensity_settings$rows)) testL <- apollo_insertRows(testL, normalDensity_settings$rows, 1)
    if(all(testL==0)) stop("\nCALCULATION ISSUE - All observations have zero probability at starting value for model component \"", normalDensity_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("\nSome observations have zero probability at starting value for model component \"", normalDensity_settings$componentName,"\""), type="i")
    return(invisible(testL))
  }

  # ####################################### #
  #### functionality="zero_LL/shares_LL" ####
  # ####################################### #
  
  if(functionality %in% c("zero_LL","shares_LL")){
    P <- rep(NA, normalDensity_settings$nObs)
    if(any(!normalDensity_settings$rows)) P <- apollo_insertRows(P, normalDensity_settings$rows, 1)
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction"){
    ans <- normalDensity_settings$mu + normalDensity_settings$xNormal
    if(any(!normalDensity_settings$rows)) ans <- apollo_insertRows(ans, normalDensity_settings$rows, NA)
    return(ans)
  } 

  # ################################################################# #
  #### functionality="estimate/conditionals/raw/output/components" ####
  # ################################################################# #

  if(functionality %in% c("estimate", "conditionals", "raw", "output", "components")){
    ans <- stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                        normalDensity_settings$mu, normalDensity_settings$sigma)
    if(any(!normalDensity_settings$rows)){
      if(functionality=="raw"){
        ans <- apollo_insertRows(ans, normalDensity_settings$rows, NA)
      } else {
        ans <- apollo_insertRows(ans, normalDensity_settings$rows, 1)
      } 
    }
    return(ans)
  }
  
  # ################################ #
  #### functionality="gradient"   ####
  # ################################ #
  
  if(functionality=="gradient"){
    if(!normalDensity_settings$gradient) stop("INTERNAL ISSUE - Analytical gradient could not be calculated for ", 
                                              normalDensity_settings$componentName, 
                                              ". Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_normalDensity could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_normalDensity could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate likelihood
    L <- stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                      normalDensity_settings$mu, normalDensity_settings$sigma)
    # Calculate dL/db
    G <- list()
    K <- length(normalDensity_settings$dX) # number of parameters
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    eHatSig <- normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal - normalDensity_settings$mu
    eHatSig <- eHatSig/normalDensity_settings$sigma
    r <- all(normalDensity_settings$rows)
    for(k in 1:K){
      # Evaluate partial derivatives
      dX     <- normalDensity_settings$dX[[k]]    ; environment(dX    ) <- e; dX     <- dX()
      dMu    <- normalDensity_settings$dMu[[k]]   ; environment(dMu   ) <- e; dMu    <- dMu()
      dSigma <- normalDensity_settings$dSigma[[k]]; environment(dSigma) <- e; dSigma <- dSigma()
      if(!r){ # remove excluded rows
        dX     <- apollo_keepRows(dX    , normalDensity_settings$rows)
        dMu    <- apollo_keepRows(dMu   , normalDensity_settings$rows)
        dSigma <- apollo_keepRows(dSigma, normalDensity_settings$rows)
      }
      G[[k]] <- L/normalDensity_settings$sigma*( eHatSig*(dX + dMu) + dSigma*(eHatSig^2 - 1) )
    }; rm(dX, dMu, dSigma)
    
    # Restore rows
    if(!r){
      L <- apollo_insertRows(L, normalDensity_settings$rows, 1)
      G <- lapply(G, apollo_insertRows, r=normalDensity_settings$rows, val=0)
    }
    return(list(like=L, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(normalDensity_settings$normalDensity_diagnostics(normalDensity_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(normalDensity_settings$normalDensity_diagnostics(normalDensity_settings, apollo_inputs, data =FALSE))
    return(P)
  }

}
