#' Calculates density from a Normal distribution
#'
#' Calculates density from a Normal distribution at a specific value with a specified mean and standard deviation.
#'
#' This function estimates the linear model outcomeNormal = mu + xNormal + epsilon, where epsilon is a random error distributed Normal(0,sigma).
#' If using this function in the context of an Integrated Choice and Latent Variable (ICLV) model with continuous
#' indicators, then \code{outcomeNormal} would be the value of the indicator, \code{xNormal} would be the value of the latent variable (possibly
#' multiplied by a parameter to measure its correlation with the indicator, e.g. xNormal=lambda*LV), and \code{mu} would be
#' an additional parameter to be estimated (the mean of the indicator, which should be fixed to zero if the indicator is
#' centered around its mean beforehand).
#' @param normalDensity_settings List of arguments to the functions. It must contain the following.
#'                               \itemize{
#'                                 \item \code{outcomeNormal}: Numeric vector. Dependant variable.
#'                                 \item \code{xNormal}: Numeric vector. Single explanatory variable.
#'                                 \item \code{mu}: Numeric scalar. Intercept of the linear model.
#'                                 \item \code{sigma}: Numeric scalar. Variance of error component of linear model to be estimated.
#'                                 \item \code{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                                 \item \code{componentName}: Character. Name given to model component.
#'                               }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate": Used for model estimation.
#'                        \item "prediction": Used for model predictions.
#'                        \item "validate": Used for validating input.
#'                        \item "zero_LL": Used for calculating null likelihood.
#'                        \item "conditionals": Used for calculating conditionals.
#'                        \item "output": Used for preparing output after model estimation.
#'                        \item "raw": Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the likelihood for each observation.
#'           \item \strong{\code{"prediction"}}: Not implemented. Returns NA.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns NA.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
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
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", normalDensity_settings$componentName,
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
      tmp <- list(xNormal= normalDensity_settings$xNormal,
                  mu     = normalDensity_settings$mu,
                  sigma  = normalDensity_settings$sigma)
      normalDensity_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, tmp)
      normalDensity_settings$gradient <- !is.null(normalDensity_settings$dV)
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
    if(!apollo_control$noDiagnostics) apollo_diagnostics(normalDensity_settings, modelType, apollo_inputs)
    testL = stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                         normalDensity_settings$mu, normalDensity_settings$sigma)
    if(any(!normalDensity_settings$rows)) testL <- apollo_insertRows(testL, normalDensity_settings$rows, 1)
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"", normalDensity_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("\nSome observations have zero probability at starting value for model component \"", normalDensity_settings$componentName,"\""))
    return(invisible(testL))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL") return(NA)

  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction") return(NA)

  # ###################################################### #
  #### functionality="estimate/conditionals/raw/output" ####
  # ###################################################### #

  if(functionality %in% c("estimate", "conditionals", "raw", "output")){
    ans <- stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                        normalDensity_settings$mu, normalDensity_settings$sigma)
    if(any(!normalDensity_settings$rows)) ans <- apollo_insertRows(ans, normalDensity_settings$rows, 1)
    return(ans)
  }
  
  # ################################ #
  #### functionality="gradient"   ####
  # ################################ #
  
  if(functionality=="gradient"){
    if(!normalDensity_settings$gradient) stop("Analytical gradient could not be calculated for ", 
                                              normalDensity_settings$componentName, 
                                              ". Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("apollo_mnl could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("apollo_mnl could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate likelihood and derivatives of xNormal, mu and sigma
    L <- stats::dnorm(normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal,
                      normalDensity_settings$mu, normalDensity_settings$sigma)
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs)), hash=TRUE)
    for(i in 1:length(normalDensity_settings$dV)) environment(normalDensity_settings$dV[[i]]) <- e
    dV <- lapply(normalDensity_settings$dV, function(dv) dv())
    if(!all(normalDensity_settings$rows)) for(i in 1:length(dV)) dV[[i]] <- lapply(dV[[i]], 
                                                                                   apollo_keepRows, 
                                                                                   normalDensity_settings$rows)
    # Calculate gradient
    G  <- normalDensity_settings$outcomeNormal - normalDensity_settings$xNormal - normalDensity_settings$mu
    G2 <- G^2
    G <- mapply(function(dv, dm, ds) G*(dv + dm) + ds/normalDensity_settings$sigma*(G2 - normalDensity_settings$sigma^2), 
                dV$xNormal, dV$mu, dV$sigma, SIMPLIFY=FALSE)
    G <- lapply(G, function(g) L/normalDensity_settings$sigma^2*g)
    
    # Restore rows
    if(!all(normalDensity_settings$rows)){
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
    P$data  <- capture.output(apollo_diagnostics(normalDensity_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(normalDensity_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }

}
