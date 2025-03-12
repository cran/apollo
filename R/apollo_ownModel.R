#' Calculates own model probabilities
#'
#' Receives functions or expressions for each functionality so that a 
#' user-defined model can interface with Apollo.
#' 
#' @param ownModel_settings List of arguments. Only likelihood is mandatory.
#'                          \itemize{
#'                           \item \strong{\code{gradient}}: Function or
#'                                 expression used to calculate the gradient
#'                                 of the likelihood. If not provided, Apollo 
#'                                 will attempt to calculate it automatically.
#'                           \item \strong{\code{likelihood}}: Function or
#'                                 expression used to calculate the likelihood
#'                                 of the model. Should evaluate to a vector, 
#'                                 matrix, or 3-dimensional array.
#'                           \item \strong{\code{prediction}}: Function or
#'                                 expression used to calculate the prediction
#'                                 of the model. Should evaluate to a vector, 
#'                                 matrix, or 3-dimensional array.
#'                           \item \strong{\code{report}}: List of functions or
#'                                 expressions used to produce a text report 
#'                                 summarising the input and parameter 
#'                                 estimates of the model. Should contain two 
#'                                 elements: "data" (with a summary of the 
#'                                 input data), and "param" (with a summary of 
#'                                 the estimated parameters).
#'                           \item \strong{\code{shares_LL}}: Function or
#'                                 expression used to calculate the likelihood
#'                                 of the constants-only model.
#'                           \item \strong{\code{zero_LL}}: Function or
#'                                 expression used to calculate the likelihood
#'                                 of the base model (e.g. equiprobable model).
#'                          }
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
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{mnl_settings}.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Choice overview
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @export
#' @importFrom utils capture.output
apollo_ownModel <- function(ownModel_settings, functionality){
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE, silent=FALSE, analyticGrad=TRUE)) ))
  
  ### Set or extract componentName
  modelType   = "ownModel"
  if(is.null(ownModel_settings[["componentName"]])){
    ownModel_settings[["componentName"]] = ifelse(!is.null(ownModel_settings[['componentName2']]),
                                             ownModel_settings[['componentName2']], modelType)
    test <- functionality=="validate" && ownModel_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', ownModel_settings[["componentName"]], '" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, ownModel_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", ownModel_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  if( !is.null(apollo_inputs[[paste0(ownModel_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load ownModel_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(ownModel_settings$componentName, "_settings")]]
    # The only possible pre-processing for ownModel is the gradient and 
    # modelType, so those are the only thing taken from tmp, if it exists.
    if(!is.null(tmp$modelType)) ownModel_settings$modelType <- tmp$modelType
    if(!is.null(tmp$gradient )) ownModel_settings$gradient  <- tmp$gradient
    if(!is.null(tmp$dLike    )) ownModel_settings$dLike     <- tmp$dLike
    if(!is.null(tmp$nObs     )) ownModel_settings$nObs      <- tmp$nObs
    if(!is.null(tmp$rows     )) ownModel_settings$rows      <- tmp$rows
    if(!is.null(tmp$componentName)) ownModel_settings$componentName <- tmp$componentName
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    test <- is.null(names(ownModel_settings))
    if(test) stop("SYNTAX ISSUE - All elements inside the inputs lists for ", 
                  "model components must be named, e.g. \"ownModel_settings", 
                  "=list(likelihood=..., rows=...)\".")
    test <- is.null(ownModel_settings[["componentName"]])
    if(test) stop("SYNTAX ISSUE - The settings of at least one model component",
                  " is missing the mandatory \"componentName\" object.")
    # validate functionality
    test <- functionality %in% c("estimate","prediction","validate","zero_LL",
                                 "shares_LL","conditionals","output","raw",
                                 "preprocess","components","gradient","report")
    if(!test) stop("SYNTAX ISSUE - Non-permissable value of ",
                   "\"functionality\" for model component \"",
                   ownModel_settings$componentName,"\".")
    # Validate mandatory inputs
    for(i in c("likelihood")){
      test <- i %in% names(ownModel_settings)
      if(!test) stop("SYNTAX ISSUE - The inputs list for model component \"", 
                     ownModel_settings$componentName, "\" needs to include ", 
                     "an object called \"", i,"\"!")
    }
    # Validate optional inputs
    for(i in c("prediction", "zero_LL", "shares_LL", "gradient", "report")){
      txt1 <- "SYNTAX ISSUE - Argument \""
      txt3 <- paste0("\" in model component ", ownModel_settings$componentName, 
                     " must be ")
      if(i %in% names(ownModel_settings)){
        if(i %in% c("prediction", "zero_LL", "shares_LL")){
          test <- is.function(ownModel_settings[[i]])
          test <- test || is.numeric(ownModel_settings[[i]])
          test <- test || is.logical(ownModel_settings[[i]])
          if(!test) stop(txt1, i, txt3, "a function or an expression ", 
                         "evaluating to a numeric value.")
        }
        if(i=="gradient"){
          K <- tryCatch(length(get("apollo_beta", envir=parent.frame(), 
                                   inherits=FALSE)), 
                        error=function(e) 0)
          K <- K - length(apollo_inputs$apollo_fixed)
          test <- is.list(ownModel_settings[[i]]) && 
            all(is.function(ownModel_settings[[i]]))
          if(K>0) test <- test && length(ownModel_settings[[i]])==K
          if(!test) stop(txt1, i, txt3, "a list of functions, one per ", 
                         "parameter to be estimated.")
          ownModel_settings$dLike    <- ownModel_settings$gradient
          ownModel_settings$gradient <- TRUE
        }
        if(i=="report"){
          test <- is.list(ownModel_settings[[i]])
          test <- test && !is.null(names(ownModel_settings[[i]]))
          test <- test && all(c("data", "param") %in% names(ownModel_settings[[i]]))
          test <- test && all(sapply(ownModel_settings[[i]],
                                     function(x) is.function(x) | is.character(x)))
          if(!test) stop(txt1, i, txt3, "a list of two functions or expressions ", 
                         "called \"data\" and \"param\", each producing or ", 
                         "containing a description of the model data and ", 
                         "estimated parameters, respectively.")
        }
      }
    }
    # Validate rows, and set nObs and modelType
    ownModel_settings$rows <- aux_validateRows(ownModel_settings$rows, 
                                               ownModel_settings$componentName, 
                                               apollo_inputs)
    ownModel_settings$nObs <- sum(ownModel_settings$rows)
    ownModel_settings$modelType <- modelType
    # Construct necessary input for gradient
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && is.null(ownModel_settings$dLike)
    test <- test && (is.null(ownModel_settings$gradient) || !ownModel_settings$gradient)
    test <- test && is.function(ownModel_settings$likelihood)
    test <- test && apollo_inputs$apollo_control$analyticGrad
    if(is.null(ownModel_settings$gradient)) ownModel_settings$gradient <- FALSE
    if(test){
      ownModel_settings$dLike    <- apollo_dVdB(apollo_beta, apollo_inputs, 
                                                 list(ownModel_settings$likelihood))
      if(is.list(ownModel_settings$dLike)){
        ownModel_settings$dLike <- lapply(ownModel_settings$dLike, "[[", 1)
      }
      ownModel_settings$gradient <- !is.null(ownModel_settings$dLike)
    }; rm(test)
    # Return ownModel_settings if pre-processing
    if(functionality=="preprocess"){
      # Return only things that remain the same across iterations
      return(ownModel_settings[c("modelType", "gradient", "dLike", "nObs", 
                                 "rows", "componentName")])
    }
  }
  
  # ##################################################### #
  #### Transform likelihood into numeric and drop rows ####
  # ##################################################### #
  
  ### Execute likelihood (makes sure we are now working with 
    # vectors/matrices/arrays and not functions)
  test <- is.function(ownModel_settings$likelihood)
  if(test) ownModel_settings$likelihood = ownModel_settings$likelihood()
  test <- is.matrix(ownModel_settings$likelihood)
  test <- test && ncol(ownModel_settings$likelihood)==1
  if(test) ownModel_settings$likelihood = as.vector(ownModel_settings$likelihood)
  # Drops rows if necessary
  if(!all(ownModel_settings$rows)){
    ownModel_settings$likelihood <- apollo_keepRows(ownModel_settings$likelihood,
                                                    ownModel_settings$rows)
  }

  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation){
      if(ownModel_settings$nObs==0) stop("INPUT ISSUE - No data for model component \"",
                                         ownModel_settings$componentName,"\"")
    }
    
    testL = ownModel_settings$likelihood
    if(!all(ownModel_settings$rows)) testL <- apollo_insertRows(testL, ownModel_settings$rows, 1)
    if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",ownModel_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",ownModel_settings$componentName,"\"", sep=""))
    return(invisible(testL))
  }
  
  # ###################################################### #
  #### functionality="prediction/zero_LL/shares_LL/raw" ####
  # ###################################################### #
  if(functionality=="raw") functionality <- "prediction"
  
  if(functionality %in% c("zero_LL", "shares_LL", "prediction")){
    if(functionality %in% names(ownModel_settings)){
      P <- ownModel_settings[[functionality]]
      if(is.function(P)) P <- P()
      if(!all(ownModel_settings$rows)){
        P <- apollo_keepRows(P, ownModel_settings$rows)
        P <- apollo_insertRows(P, ownModel_settings$rows, 1)
      }
      return(P)
    } else return(NULL)
  }
  
  # ############################################################# #
  #### functionality="estimate/prediction/conditionals/output" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","conditionals", "components", "output")){
    P <- ownModel_settings$likelihood
    if(!all(ownModel_settings$rows)) P <- apollo_insertRows(P, ownModel_settings$rows, 1)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    test <- !is.null(ownModel_settings$dLike) && is.list(ownModel_settings$dLike)
    test <- test && all(sapply(ownModel_settings$dLike, is.function))
    if(!test) stop("CALCULATION ISSUE - No analytical gradient available")
    
    # Evaliate gradient
    G <- lapply(ownModel_settings$dLike, function(f) f())
    P <- ownModel_settings$likelihood
    # Remove and re-insert rows if necessary
    if(!all(ownModel_settings$rows)){
      G <- lapply(G, apollo_keepRows  , r=ownModel_settings$rows)
      G <- lapply(G, apollo_insertRows, r=ownModel_settings$rows, val=1)
      P <- apollo_insertRows(P, r = ownModel_settings$rows, val = 1)
    }
    return(list(like=P, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    if(functionality %in% names(ownModel_settings)){
      P <- ownModel_settings[[functionality]]
      P <- lapply(P, function(p) if(is.function(p)) p() else p)
      return(P)
    } else return(NULL)
  }
  
}
