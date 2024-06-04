#' Creates hessian function.
#'
#' Creates hessian function from the likelihood function apollo_probabilities provided by the user. Returns NULL if 
#' the creation of gradient function fails.
#'
#' Internal use only. Called by \code{apollo_estimate} before estimation.
#' The returned function can be single-threaded or multi-threaded based on the model options.
#' @param apollo_beta Named numeric vector. Names and values for (all) parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not 
#'                     change during estimation.
#' @param apollo_logLike Function to calculate the log-likelihood of the model, as created by \link{apollo_makeLogLike}
#'                       If provided, the value of the analytical gradient will be compared to the value of the
#'                       numerical gradient as calculated using apollo_logLike and the numDeriv package.
#'                       If the difference between the two is bigger than 1% for any dimension, it will be assumed
#'                       that the analytical gradient is wrong and NULL will be returned.
#' @return apollo_hessian function. It receives a single argument called \code{b}, which are the _variable_ 
#'         parameters (i.e. must not include fixed parameters).
#' @export
apollo_makeHessian <- function(apollo_beta, apollo_fixed, apollo_logLike){
  
  # # # # # # # # # # # # #
  #### Fetch variables ####
  # # # # # # # # # # # # #
  
  # apollo_probabilities, apollo_inputs, singleCore, debug, and silent 
  if(!is.null(environment(apollo_logLike)[['cl']])) cl <- environment(apollo_logLike)[['cl']] else cl <- NA
  singleCore <- !is.list(cl) || (length(cl)==1 && is.na(cl))
  if(singleCore){
    apollo_probabilities <- environment(apollo_logLike)[['apollo_probabilities']]
    apollo_inputs <- environment(apollo_logLike)[['apollo_inputs']]
  } else {
    apollo_probabilities <- parallel::clusterEvalQ(cl, apollo_probabilities)[[1]]
    apollo_inputs <- parallel::clusterEvalQ(cl, apollo_inputs[-which(names(apollo_inputs) %in% c('database', 'draws'))])[[1]]
  }
  if(!is.null(apollo_inputs$apollo_control$debug)) debug <- apollo_inputs$apollo_control$debug else debug <- FALSE
  if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
  
  
  # # # # # # # # #
  #### Checks ####
  # # # # # # # # #
  
  ### Check that no models without analytical hessian are used in apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    txt <- c("apollo_fmnl|apollo_op|apollo_normalDensity|apollo_dft|apollo_mdcev|apollo_el|apollo_nl|apollo_cnl|apollo_mdcnev|apollo_op|apollo_emdc1|apollo_emdc2|apollo_emdc|apollo_fnl")
    tmp <- grep(txt, tmp)
    if(length(tmp)>0){
      if(debug) apollo_print("Analytic gradient cannot be built because models with undefined gradient are used inside apollo_probabilities.")
      return(NULL)
    }; rm(txt, tmp)
  }
  
  ### Turn off analytic hessian if using inter-intra, unless manually set to TRUE
  if(apollo_inputs$apollo_control$mixing && is.list(apollo_inputs$apollo_draws)){
    test <- apollo_inputs$apollo_draws$interNDraws>1
    test <- test & apollo_inputs$apollo_draws$intraNDraws>1
    test <- test & !is.null(apollo_inputs$apollo_control$analyticGrad_manualSet) 
    test <- test && !apollo_inputs$apollo_control$analyticGrad_manualSet
    if(test & debug) apollo_print("By default, analytic hessian will not be used for your model as it combines inter & intra draws (to avoid excesive memory usage). If you wish to use analytic gradients, please set analyticHessian = TRUE in apollo_control.")
    if(test) return(NULL)
  }
  
  ### Check that second derivatives are available for all components
  if(!singleCore){ # multi-core
    d2VAvail <- parallel::clusterEvalQ(cl, {
      compNames <- grep("_settings$", names(apollo_inputs), value=TRUE)
      d2VAvail   <- c()
      for(i in compNames) d2VAvail <- c(d2VAvail, apollo_inputs[[i]]$hessian)
      compNames <- substr(compNames, 1, nchar(compNames)-nchar("_settings"))
      if(length(compNames)>0) setNames(d2VAvail, compNames) else d2VAvail
    })[[1]]
  } else { # single-core
    compNames <- grep("_settings$", names(apollo_inputs), value=TRUE)
    d2VAvail   <- c()
    for(i in compNames) d2VAvail <- c(d2VAvail, apollo_inputs[[i]]$hessian)
    compNames <- substr(compNames, 1, nchar(compNames)-nchar("_settings"))
    if(length(compNames)>0) d2VAvail <- setNames(d2VAvail, compNames)
    rm(compNames, i)
  }
  if( length(d2VAvail)==0 || any(!d2VAvail)){
    txt <- paste('Apollo was not able to compute analytical hessians for your',
                 'model. This could be because you are using model components', 
                 'for which analytical hessians are not yet implemented, or',
                 'because you coded your own model functions. If however you',
                 'only used apollo_mnl, then there could be another issue.',
                 'You might want to ask for help in the Apollo forum',
                 '(http://www.apollochoicemodelling.com/forum) on how to', 
                 'solve this issue. If you do, please post your code and',
                 'data (if not confidential).')
    if(!silent) apollo_print(txt,  pause=0, type="i")
    return(NULL)
  }
  if( !all(d2VAvail) ) return(NULL)
  
  
  # # # # # # # # # # # # # # # # #
  #### Build Hessian function ####
  # # # # # # # # # # # # # # # # #
  
  ### Split apollo_beta in fixed and variable parts
  bFix <- apollo_beta[apollo_fixed]
  bOrd <- names(apollo_beta)
  
  ### Construct gradient function
  if(singleCore){ # Single core
    hessian <- function(b){
      b <- c(b, bFix)[bOrd]
      H <- apollo_probabilities(b, apollo_inputs, functionality="hessian")
      return( H )
    }
    environment(hessian) <- environment(apollo_logLike)
  } else { # Multi-core
    hessian <- function(b){
      b <- c(b, bFix)[bOrd]
      parallel::clusterExport(cl, "b", envir=environment())
      H <- parallel::clusterEvalQ(cl=cl, apollo_probabilities(b, apollo_inputs, functionality="hessian"))
      H <- do.call("+", H)
      return( H )
    }
    environment(hessian) <- new.env(parent=baseenv())
    assign("cl", cl, envir=environment(hessian))
  }
  
  ### Copy elements to function environment
  assign("bFix", bFix, envir=environment(hessian))
  assign("bOrd", bOrd, envir=environment(hessian))
  assign("silent"     ,     silent, envir=environment(hessian))
  assign("debug"      ,      debug, envir=environment(hessian))
  assign("singleCore" , singleCore, envir=environment(hessian))
  
  
  return(hessian)
}