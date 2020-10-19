#' Creates gradient function.
#'
#' Creates gradient function from the likelihood function apollo_probabilities provided by the user. Returns NULL if 
#' the creation of gradient function fails.
#'
#' Internal use only. Called by \code{apollo_estimate} before estimation.
#' The returned function can be single-threaded or multi-threaded based on the model options.
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not 
#'                     change during estimation.
#' @param apollo_logLike Function to calculate the loglikelihood of the model, as created by \link{apollo_makeLogLike}
#'                       If provided, the value of the analytical gradient will be compared to the value of the
#'                       numerical gradient as calculated using apollo_logLike and the numDeriv package.
#'                       If the difference between the two is bigger than 1% for any dimension, it will be assumed
#'                       that the analytical gradient is wrong and NULL will be returned.
#' @param validateGrad Logical. If TRUE, it compares the value of the analytical gradient evaluated at apollo_beta
#'                     against the numeric gradient (using numDeriv) at the same value. If the difference is bigger
#'                     than 1% for any dimension, it will assume that the analytical gradient is wrong and will
#'                     return NULL.
#' @return apollo_gradient function. It receives the following arguments
#'         \itemize{
#'           \item \code{b} Numeric vector of _variable_ parameters (i.e. must not include fixed parameters).
#'           \item \code{countIter} Not used. Included only to mirror inputs of apollo_logLike.
#'           \item \code{writeIter} Not used. Included only to mirror inputs of apollo_logLike.
#'           \item \code{sumLL} Not used. Included only to mirror inputs of apollo_logLike.
#'           \item \code{getNIter} Not used. Included only to mirror inputs of apollo_logLike.
#'         }
#'         If the creation of the gradient function fails, then it returns NULL.
#' @export
apollo_makeGrad <- function(apollo_beta, apollo_fixed, apollo_logLike, validateGrad=FALSE){
  
  ### Extract variables
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
  
  # # # # # #  # #
  #### Checks ####
  # # # # # #  # #
  
  ### Check that no models without analytical gradient are used in apollo_probabilities
  if(is.function(apollo_probabilities)){
    tmp <- as.character(body(apollo_probabilities))
    txt <- c("apollo_lc|apollo_dft|apollo_mdcev|apollo_el|apollo_nl|apollo_cnl|apollo_mdcnev")
    tmp <- grep(txt, tmp)
    if(length(tmp)>0){
      if(debug) apollo_print("Analytic gradient cannot be built because models with undefined gradient are used inside apollo_probabilities.")
      return(NULL)
    }
    rm(txt, tmp)
  }
  
  ### Turn off analytic gradient if using inter-intra
  if(apollo_inputs$apollo_control$mixing && is.list(apollo_inputs$apollo_draws)){
    test <- apollo_inputs$apollo_draws$interNDraws>1
    test <- test & apollo_inputs$apollo_draws$intraNDraws>1
    if(test & debug) apollo_print("Analytic gradients will not be built because inter & intra draws are being used (to avoid excesive memory usage).")
    if(test) return(NULL)
  }
  
  ### Check that gradients are available for all components
  if(!singleCore){ # multi-core
    dVAvail <- parallel::clusterEvalQ(cl, {
      compNames <- grep("_settings$", ls(apollo_inputs), value=TRUE)
      dVAvail   <- c()
      for(i in compNames) dVAvail <- c(dVAvail, apollo_inputs[[i]]$gradient)
      compNames <- substr(compNames, 1, nchar(compNames)-nchar("_settings"))
      if(length(compNames)>0) setNames(dVAvail, compNames) else dVAvail
    })[[1]]
  } else { # single-core
    compNames <- grep("_settings$", ls(apollo_inputs), value=TRUE)
    dVAvail   <- c()
    for(i in compNames) dVAvail <- c(dVAvail, apollo_inputs[[i]]$gradient)
    compNames <- substr(compNames, 1, nchar(compNames)-nchar("_settings"))
    if(length(compNames)>0) dVAvail <- setNames(dVAvail, compNames)
    rm(compNames, i)
  }
  compList <- apollo_inputs$apolloLog$listOfNames
  if( length(dVAvail)==0 || !all(compList %in% names(dVAvail)) ){
    txt <- paste0("Some model components cannot be pre-processed. ",
                  "Numeric gradients will be used")
    if(debug) apollo_print(txt)
    return(NULL)
  } 
  if( !all(dVAvail) ) return(NULL)
  rm(compList)
  
  
  # # # # # # # # # # # # # # # # #
  #### Build gradient function ####
  # # # # # # # # # # # # # # # # #
  
  ### Split apollo_beta in fixed and variable parts
  bFix <- apollo_beta[apollo_fixed]
  bOrd <- names(apollo_beta)
  
  ### Construct gradient function
  if(singleCore){ # Single core
    grad <- function(b, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      b <- c(b, bFix)[bOrd]
      G <- apollo_probabilities(b, apollo_inputs, functionality="gradient")
      return( G )
    }
    #environment(grad) <- new.env(parent=baseenv()) # parent used to be parent.frame()
    #assign('apollo_inputs', apollo_inputs, envir=environment(grad))
    #assign('apollo_probabilities', apollo_probabilities, envir=environment(grad))
    environment(grad) <- environment(apollo_logLike)
  } else { # Multi-core
    grad <- function(b, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      b <- c(b, bFix)[bOrd]
      parallel::clusterExport(cl, "b", envir=environment())
      G <- parallel::clusterEvalQ(cl=cl, apollo_probabilities(b, apollo_inputs, functionality="gradient"))
      G <- do.call(rbind, G)
      return( G )
    }
    environment(grad) <- new.env(parent=baseenv())
    assign("cl", cl, envir=environment(grad))
  }
  
  ### Copy elements to function environment
  assign("bFix", bFix, envir=environment(grad))
  assign("bOrd", bOrd, envir=environment(grad))
  assign("silent"      ,     silent, envir=environment(grad))
  assign("debug"       ,      debug, envir=environment(grad))
  assign("singleCore" ,  singleCore, envir=environment(grad))
  
  
  # # # # # # # # # # # # # # # # #
  #### Test gradient function  ####
  # # # # # # # # # # # # # # # # #
  
  if(!is.null(grad) && (validateGrad || !is.null(apollo_logLike))){
    if(debug) apollo_print(c("\n", "Validating gradient function..."))
    b0_disturbance <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
    
    # Calculate analytical gradient
    gradAn <- tryCatch(grad(b0_disturbance), error=function(e) return(NULL))
    if(is.null(gradAn) || anyNA(gradAn)){ # If analytical gradient failed
      if(debug) apollo_print("Analytic gradient calculation failed. Defaulting to numeric one.")
      if(debug && !is.null(gradAn) && anyNA(gradAn)){
        colsNA <- which(apply(gradAn, MARGIN=2, anyNA))
        if(apollo_inputs$apollo_control$panelData){
          apollo_print("Parameters and individuals with NA in the analytic gradient:")
          if(singleCore) ind <- unique(apollo_inputs$database[, apollo_inputs$apollo_control$indivID])
          if(!singleCore){
            ind <- parallel::clusterEvalQ(cl, apollo_inputs$database[, apollo_inputs$apollo_control$indivID])
            ind <- unique(unlist(ind))
          } 
          for(i in colsNA){
            tmp <- ind[which(is.na(gradAn[,i]))]
            if(length(tmp)>30) tmp <- c(as.character(tmp[1:30]), paste0('... (', length(tmp), ' indivs in total)'))
            apollo_print(paste0(names(b0_disturbance)[i], ':', paste0(tmp, collapse=", ")))
          }
        } else {
          apollo_print("Parameters and rows in database with NA in the analytic gradient:")
          for(i in colsNA){
            tmp <- which(is.na(gradAn[,i]))
            if(length(tmp)>50) tmp <- c(as.character(tmp[1:50]), paste0('... (', length(tmp), ' rows in total)'))
            apollo_print(paste0(names(b0_disturbance)[i], ':', paste0(tmp, collapse=", ")))
          }
        }
      }
      return(NULL)
    } else gradAn <- colSums(gradAn)
    
    if(!is.null(gradAn) && all(is.finite(gradAn)) && is.function(apollo_logLike)){
      # Calculate numerical gradient for testing
      gradNum <- tryCatch(numDeriv::grad(apollo_logLike, b0_disturbance, sumLL=TRUE), error=function(e) return(NULL))
      if(!is.null(gradNum)) names(gradNum) <- names(b0_disturbance)
      if(is.null(gradNum) && debug) apollo_print("Numeric gradient calculation failed.")
      
      # If analytical and numeric gradient calculation was succesful
      test <- !is.null(gradNum) && all(is.finite(gradNum))
      dif  <- abs((abs(gradNum)-abs(gradAn))/abs(gradNum))
      test <- test && any(dif[is.finite(dif)]>0.01) # diff should be <1% for all elements
      if(test){
        if(debug) { apollo_print("Gradient pre optimisation : "); print(gradNum) }
        if(debug) { apollo_print("Gradient post optimisation: "); print(gradAn ) }
        if(!silent) apollo_print("Analytical gradient is different to numerical one. Numerical gradients will be used.")
        return(NULL)
      }; rm(test)
    }
    
  }; rm(b0_disturbance)
  
  
  return(grad)
}