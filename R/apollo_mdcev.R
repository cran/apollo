#' Calculates MDCEV likelihoods
#'
#' Calculates the likelihoods of a Multiple Discrete Continuous Extreme Value (MDCEV) model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' @param mdcev_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                       \itemize{
#'                         \item \strong{\code{alpha}}: Named list. Alpha parameters for each alternative, including for any outside good. As many elements as alternatives.
#'                         \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'utilities'.
#'                         \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                         \item \strong{\code{budget}}: Numeric vector. Budget for each observation.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                         \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{cost}}: Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{gamma}}: Named list. Gamma parameters for each alternative, excluding any outside good. As many elements as inside good alternatives.
#'                         \item \strong{\code{nRep}}: Numeric scalar. Number of simulations of the whole dataset used for forecasting. The forecast is the average of these simulations. Default is 100.
#'                         \item \strong{\code{outside}}: Character. Optional name of the outside good.
#'                         \item \strong{\code{rawPrediction}}: Logical scalar. TRUE for prediction to be returned at the draw level (a 3-dim array). FALSE for prediction to be returned averaged across draws. Default is FALSE.
#'                         \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                         \item \strong{\code{sigma}}: Numeric scalar. Scale parameter of the model extreme value type I error.
#'                         \item \strong{\code{utilities}}: Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#'                       }
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
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the observed consumption for each observation.
#'           \item \strong{\code{"gradient"}}: Not implemented
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: A matrix with one row per observation, and columns indicating means and s.d. of continuous and discrete predicted consumptions.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{mdcev_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"report"}}: Dependent variable overview.
#'           \item \strong{\code{"shares_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'         }
#' @importFrom utils capture.output
#' @export
apollo_mdcev <- function(mdcev_settings,functionality){
  ### Set or extract componentName
  modelType   = "MDCEV"
  if(is.null(mdcev_settings[["componentName"]])){
    mdcev_settings[["componentName"]] = ifelse(!is.null(mdcev_settings[['componentName2']]),
                                               mdcev_settings[['componentName2']], modelType)
    test <- functionality=="validate" && mdcev_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 mdcev_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, mdcev_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", mdcev_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used
  if(!is.null(mdcev_settings[["utilities"]])) names(mdcev_settings)[which(names(mdcev_settings)=="utilities")]="V"
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(mdcev_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    
    # Load mdcev_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(mdcev_settings$componentName, "_settings")]]
    # If there is no V inside the loaded mdcev_settings, restore the one received as argument
    if(is.null(tmp$V    )) tmp$V     <- mdcev_settings$V    
    if(is.null(tmp$alpha)) tmp$alpha <- mdcev_settings$alpha
    if(is.null(tmp$gamma)) tmp$gamma <- mdcev_settings$gamma
    if(is.null(tmp$sigma)) tmp$sigma <- mdcev_settings$sigma
    mdcev_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    mdcev_settings <- apollo_preprocess(mdcev_settings, modelType, functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation available for OL components.")
    # Using R likelihood
    if(mdcev_settings$hasOutside) mdcev_settings$probs_MDCEV <- function(inputs){
      # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues alpha
      inputs$V     <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V    , inputs$avail, SIMPLIFY=FALSE)
      inputs$alpha <- mapply(function(l,a) apollo_setRows(l, !a, 0), inputs$alpha, inputs$avail, SIMPLIFY=FALSE)
      inputs$gamma <- mapply(function(g,a) apollo_setRows(g, !a, 0), inputs$gamma, inputs$avail, SIMPLIFY=FALSE)
      # Compute V
      inputs$V[[1]]=(inputs$alpha[[1]]-1)*log(inputs$continuousChoice[[1]])
      for(j in 2:inputs$nAlt){
        if(inputs$minX){
          tmp <- inputs$continuousChoice[[j]]-(inputs$continuousChoice[[j]]>=inputs$minConsumption[[j]])*inputs$minConsumption[[j]]
          inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]]*((inputs$alpha[[j]]-1)*log((tmp/inputs$gamma[[j]]) + 1) - log(inputs$cost[[j]]))
        } else {
          inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]]*((inputs$alpha[[j]]-1)*log((inputs$continuousChoice[[j]]/inputs$gamma[[j]]) + 1) - log(inputs$cost[[j]]))
        }
      }
      # First term
      term1=(1-inputs$totalChosen)*log(inputs$sigma)
      # Second term
      logfi=list()
      for(j in 2:inputs$nAlt){
        if(inputs$minX){
          tmp <- inputs$continuousChoice[[j]]-(inputs$continuousChoice[[j]]>=inputs$minConsumption[[j]])*inputs$minConsumption[[j]]
          logfi[[j-1]]=inputs$avail[[j]]*( log(1-inputs$alpha[[j]]) - log(tmp + inputs$gamma[[j]]) )
        } else {
          logfi[[j-1]]=inputs$avail[[j]]*(log(1-inputs$alpha[[j]]) - log(inputs$continuousChoice[[j]] + inputs$gamma[[j]]))
        }
      }
      term2 = log(1-inputs$alpha[[1]])-log(inputs$continuousChoice[[1]])
      for(j in 2:inputs$nAlt) term2 = term2 + inputs$avail[[j]]*(logfi[[j-1]]*inputs$discrete_choice[[j]])
      # Third term
      term3 = inputs$continuousChoice[[1]]/(1-inputs$alpha[[1]])
      for(j in 2:inputs$nAlt) term3 = term3 + inputs$avail[[j]]*(inputs$cost[[j]]/exp(logfi[[j-1]])*inputs$discrete_choice[[j]])
      term3 = log(term3)
      # Fourth term
      term4_1 = inputs$V[[1]]/inputs$sigma
      term4_2 = exp(inputs$V[[1]]/inputs$sigma)
      for(j in 2:inputs$nAlt){
        term4_1 = term4_1 + inputs$avail[[j]]*(inputs$V[[j]]/inputs$sigma * inputs$discrete_choice[[j]])
        term4_2 = term4_2 + inputs$avail[[j]]*exp(inputs$V[[j]]/inputs$sigma)
      }
      term4_2 = inputs$totalChosen * log(term4_2) 
      term4 =  term4_1 - term4_2 
      rm(term4_1, term4_2)
      # Fifth term: log of factorial
      term5 = lfactorial(inputs$totalChosen-1)
      # probability is simply the exp of the sum of the logs of the individual terms
      P = exp(term1 + term2 + term3 + term4 + term5)
      rm(term1, term2, term3, term4, term5)
      # If an aunavailable alternative was chosen, likelihood is zero
      if(any(inputs$chosenUnavail)) P <- apollo_setRows(P, inputs$chosenUnavail, 0)
      return(P)
    }
    if(!mdcev_settings$hasOutside) mdcev_settings$probs_MDCEV <- function(inputs){
      # Compute V
      for(j in 1:inputs$nAlt){
        if(inputs$minX){
          tmp <- inputs$continuousChoice[[j]]-(inputs$continuousChoice[[j]]>=inputs$minConsumption[[j]])*inputs$minConsumption[[j]]
          inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]]*((inputs$alpha[[j]]-1)*log((tmp/inputs$gamma[[j]]) + 1) - log(inputs$cost[[j]]))
        } else {
          inputs$V[[j]] = inputs$V[[j]] + inputs$avail[[j]]*((inputs$alpha[[j]]-1)*log((inputs$continuousChoice[[j]]/inputs$gamma[[j]]) + 1) - log(inputs$cost[[j]]))
        }
      }
      # first term
      term1 = (1-inputs$totalChosen)*log(inputs$sigma)
      # compute log of fi for inside goods for use in second term, this has one column per inside good
      logfi = list()
      for(j in 1:inputs$nAlt){
        if(inputs$minX){
          tmp <- inputs$continuousChoice[[j]]-(inputs$continuousChoice[[j]]>=inputs$minConsumption[[j]])*inputs$minConsumption[[j]]
          logfi[[j]]=inputs$avail[[j]]*( log(1-inputs$alpha[[j]]) - log(tmp + inputs$gamma[[j]]) )     
        } else {
          logfi[[j]]=inputs$avail[[j]]*(log(1-inputs$alpha[[j]])-log(inputs$continuousChoice[[j]] + inputs$gamma[[j]]))  
        }
      }
      # second term
      term2 = 0  
      for(j in 1:inputs$nAlt) term2 = term2 + inputs$avail[[j]]*(logfi[[j]]*inputs$discrete_choice[[j]])
      # Third term
      term3 = 0  
      for(j in 1:inputs$nAlt) term3 = term3 + inputs$avail[[j]]*(inputs$cost[[j]]/exp(logfi[[j]])*inputs$discrete_choice[[j]])
      term3 = log(term3)
      # Fourth term
      term4_1 = 0
      term4_2 = 0
      for(j in 1:inputs$nAlt){
        term4_1 = term4_1 + inputs$avail[[j]]*(inputs$V[[j]]/inputs$sigma*inputs$discrete_choice[[j]])
        term4_2 = term4_2 + inputs$avail[[j]]*exp(inputs$V[[j]]/inputs$sigma)
      }
      term4_2 = inputs$totalChosen * log(term4_2) # log of the above to the power of M (totalChosen)
      term4 =  term4_1 - term4_2 # log of fourth term is now the difference of the two parts
      rm(term4_1, term4_2)
      # fifth term: log of factorial
      term5 = lfactorial(inputs$totalChosen-1)
      # probability is simply the exp of the sum of the logs of the individual terms
      P = exp(term1 + term2 + term3 + term4 + term5)
      rm(term1, term2, term3, term4, term5)
      # If an aunavailable alternative was chosen, likelihood is zero
      if(any(inputs$chosenUnavail)) P <- apollo_setRows(P, inputs$chosenUnavail, 0)
      return(P)
    }
    
    mdcev_settings$mdcev_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
      
        # Table describing dependent variable
        choicematrix <- matrix(0, nrow=4, ncol=inputs$nAlt, 
                               dimnames=list(c("Times available","Observations in which chosen",
                                               "Average consumption when available",
                                               "Average consumption when chosen"),
                                             inputs$altnames))
        for(a in 1:inputs$nAlt){
          choicematrix[1,a] <- ifelse(length(inputs$avail[[a]])==1 && inputs$avail[[a]]==1, 
                                      inputs$nObs, sum(inputs$avail[[a]]) )
          choicematrix[2,a] <- sum(inputs$discrete_choice[[a]])
          choicematrix[3,a] <- ifelse(choicematrix[1,a]>0, sum(inputs$continuousChoice[[a]])/choicematrix[1,a], 0)
          choicematrix[4,a] <- ifelse(choicematrix[2,a]>0, sum(inputs$continuousChoice[[a]])/choicematrix[2,a], 0)
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
            if(choicematrix[2,a]==choicematrix[1,a] && is.na(inputs$outside)) stop(paste0('SPECIFICATION ISSUE - Alternative "', inputs$altnames[a], '" is always chosen when available in model component "', inputs$componentName, '" and should be treated as the outside good.'))
            if(choicematrix[2,a]==choicematrix[1,a] && !is.na(inputs$outside) && inputs$altnames[a]!=inputs$outside) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always chosen when available in model component "', inputs$componentName, '".'), type="w")
          }
          #if(inputs$avail_set==TRUE & !apollo_inputs$silent) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'), type="i")
        }
        
      #### return ####
      return(invisible(TRUE))
    }
    
    
    # Store model type
    mdcev_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(mdcev_settings$V, is.function))
    test <- test && is.function(mdcev_settings$alpha)
    test <- test && is.function(mdcev_settings$gamma)
    test <- test && is.function(mdcev_settings$sigma)
    test <- test && apollo_inputs$apollo_control$analyticGrad
    mdcev_settings$gradient <- FALSE
    if(test){
      mdcev_settings$dV     <- apollo_dVdB(apollo_beta, apollo_inputs, mdcev_settings$V)
      mdcev_settings$dAlpha <- apollo_dVdB(apollo_beta, apollo_inputs, mdcev_settings$alpha)
      mdcev_settings$dGamma <- apollo_dVdB(apollo_beta, apollo_inputs, mdcev_settings$gamma)
      mdcev_settings$dSigma <- apollo_dVdB(apollo_beta, apollo_inputs, list(dSigma=mdcev_settings$sigma))[[1]]
      #mdcev_settings$gradient <- !is.null(mdcev_settings$dV)
    }; rm(test)
    
    # Return mdcev_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      mdcev_settings$V     <- NULL
      mdcev_settings$alpha <- NULL
      mdcev_settings$gamma <- NULL
      mdcev_settings$sigma <- NULL
      return(mdcev_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V, alpha, gamma and sigma (makes sure we are now working with vectors/matrices/arrays and not functions)
  testV <- any(sapply(mdcev_settings$V    , is.function))
  testA <- any(sapply(mdcev_settings$alpha, is.function))
  testG <- any(sapply(mdcev_settings$gamma, is.function))
  testS <- is.function(mdcev_settings$sigma)
  if(testV) mdcev_settings$V     = lapply(mdcev_settings$V    , function(f) if(is.function(f)) f() else f )
  if(testA) mdcev_settings$alpha = lapply(mdcev_settings$alpha, function(f) if(is.function(f)) f() else f )
  if(testG) mdcev_settings$gamma = lapply(mdcev_settings$gamma, function(f) if(is.function(f)) f() else f )
  if(testS) mdcev_settings$sigma = mdcev_settings$sigma()
  rm(testV, testA, testG, testS)
  mdcev_settings$V     <- lapply(mdcev_settings$V    , function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  mdcev_settings$alpha <- lapply(mdcev_settings$alpha, function(a) if(is.matrix(a) && ncol(a)==1) as.vector(a) else a)
  mdcev_settings$gamma <- lapply(mdcev_settings$gamma, function(g) if(is.matrix(g) && ncol(g)==1) as.vector(g) else g)
  if(is.matrix(mdcev_settings$sigma) && ncol(mdcev_settings$sigma)==1) mdcev_settings$sigma <- as.vector(mdcev_settings$sigma)
  
  ### Deal with outside
  if(mdcev_settings$hasOutside){
    # Add gamma outside, if missing
    if(is.null(mdcev_settings$gamma[[mdcev_settings$outside]])) mdcev_settings$gamma[[mdcev_settings$outside]] <- 1
    # Replace mdcev_settings$outside by "outside" in V, alpha, and gamma
    tmp <- which(names(mdcev_settings$V    )==mdcev_settings$outside); if(length(tmp)>0) names(mdcev_settings$V    )[tmp] <- "outside"
    tmp <- which(names(mdcev_settings$alpha)==mdcev_settings$outside); if(length(tmp)>0) names(mdcev_settings$alpha)[tmp] <- "outside"
    tmp <- which(names(mdcev_settings$gamma)==mdcev_settings$outside); if(length(tmp)>0) names(mdcev_settings$gamma)[tmp] <- "outside"
    rm(tmp)
  }
  ### Reorder V, alpha and gamma if necessary
  if( any(mdcev_settings$alternatives != names(mdcev_settings$V    )) ) mdcev_settings$V     <- mdcev_settings$V[mdcev_settings$alternatives]
  if( any(mdcev_settings$alternatives != names(mdcev_settings$alpha)) ) mdcev_settings$alpha <- mdcev_settings$alpha[mdcev_settings$alternatives]
  if( any(mdcev_settings$alternatives != names(mdcev_settings$gamma)) ) mdcev_settings$gamma <- mdcev_settings$gamma[mdcev_settings$alternatives]
  
  ### Reorder V and drop rows if neccesary
  if(!all(mdcev_settings$rows)){
    mdcev_settings$V     <- lapply(mdcev_settings$V    , apollo_keepRows, r=mdcev_settings$rows)
    mdcev_settings$alpha <- lapply(mdcev_settings$alpha, apollo_keepRows, r=mdcev_settings$rows)
    mdcev_settings$gamma <- lapply(mdcev_settings$gamma, apollo_keepRows, r=mdcev_settings$rows)
    mdcev_settings$sigma <- apollo_keepRows(mdcev_settings$sigma, r=mdcev_settings$rows)
    mdcev_settings$nObs  <- sum(mdcev_settings$rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    test <- !apollo_inputs$apollo_control$noValidation
    if(test) apollo_validate(mdcev_settings, modelType, functionality, apollo_inputs)
    
    test <- !apollo_inputs$apollo_control$noDiagnostics
    if(test) mdcev_settings$mdcev_diagnostics(mdcev_settings, apollo_inputs)
    
    testL = mdcev_settings$probs_MDCEV(mdcev_settings)
    if(any(!mdcev_settings$rows)) testL <- apollo_insertRows(testL, mdcev_settings$rows, 1)
    if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",mdcev_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",mdcev_settings$componentName,"\"", sep=""), type="i")
    return(invisible(testL))
  }
  
  # ####################################### #
  #### functionality="zero_LL/shares_LL" ####
  # ####################################### #
  
  if(functionality %in% c("zero_LL","shares_LL")){
    P <- rep(NA, mdcev_settings$nObs)
    if(any(!mdcev_settings$rows)) P <- apollo_insertRows(P, mdcev_settings$rows, 1)
    return(P)
  }
  
  # ############################################################# #
  #### functionality="estimate/conditionals/output/components" ####
  # ############################################################# #
  
  if(functionality %in% c("estimate","conditionals", "output", "components")){
    L <- mdcev_settings$probs_MDCEV(mdcev_settings)
    if(any(!mdcev_settings$rows)) L <- apollo_insertRows(L, mdcev_settings$rows, 1)
    return(L)
  }
  
  # #################################### #
  #### functionality="prediction/raw" ####
  # #################################### #
  
  if(functionality %in% c("prediction","raw")){
    # Change name to mdcev_settings to "s"
    s <- mdcev_settings
    rm(mdcev_settings)
    
    # Check that sigma is not random
    if(!is.vector(s$sigma)) stop("INCORRECT FUNCTION/SETTING USE - Forecasting not available for MDCEV with random sigma")
    
    # Generate draws for Gumbel error components
    if(!is.null(apollo_inputs$apollo_control$seed)) seed <- apollo_inputs$apollo_control$seed + 5 else seed <- 13 + 5
    set.seed(seed)
    tmp1 <- -log(-log(apollo_mlhs(s$nRep, s$nAlt, s$nObs)))
    epsL <- vector(mode="list", length=s$nRep)
    tmp2 <- (0:(s$nObs-1))*s$nRep
    for(r in 1:s$nRep) epsL[[r]] <- tmp1[r+tmp2,]#split(tmp1[r+tmp2,], f=rep(1:s$nAlt, each=s$nObs))
    rm(tmp1, tmp2)
    
    # Initialise output matrices
    if(s$rawPrediction) XX <- array(NA, dim=c(s$nObs, s$nAlt, s$nRep), dimnames=list(NULL, s$alternatives, NULL)) else {
      Xm <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
      Xv <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
      Mm <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
      Mv <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
      Em <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
      Ev <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    }
    X  <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    
    # Tools to deal with draws
    extractDraw <- function(b, iInter, iIntra){
      if(!is.list(b)){
        if(is.vector(b) && length(b)==1) return(rep(b, s$nObs))
        if(is.vector(b)) return(b)
        if(is.matrix(b)) return(b[,iInter])
        if(is.array(b) && length(dim(b))==3) return(b[,iInter,iIntra])
      } else {
        ans <- lapply(b, extractDraw, iInter=iInter, iIntra=iIntra)
        ans <- do.call(cbind, ans)
        return(ans)
      }
    }
    if(anyNA(apollo_inputs$apollo_draws)){
      nInter <- 0; nIntra <- 0
    } else {
      nInter <- apollo_inputs$apollo_draws$interNDraws
      nIntra <- apollo_inputs$apollo_draws$intraNDraws
    }
    if(nInter==0) nInter <- 1
    if(nIntra==0) nIntra <- 1
    step1 <- ceiling(s$nRep/10)
    step2 <- ceiling(s$nRep/2)
    
    # Expand avail, cost and budget
    avail <- sapply(s$avail, function(a) if(length(a)==1) rep(a, s$nObs) else a)
    cost  <- sapply(s$cost , function(x) if(length(x)==1) rep(x, s$nObs) else x)
    if(length(s$budget)==1) budget <- rep(s$budget, s$nObs) else budget <- s$budget
    
    ### Simulate prediction
    if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
    if(!silent) cat("0%")
    if(s$hasOutside){
      for(r in 1:s$nRep){
        X    <- 0*X
        for(iInter in 1:nInter){
          Xintra <- 0*X
          for(iIntra in 1:nIntra){
            # Extract appropiate values
            phiP <- extractDraw(s$sigma, iInter, iIntra)*epsL[[r]]
            phiP <- exp(extractDraw(s$V, iInter, iIntra) + phiP)*avail/cost
            gamma<- extractDraw(s$gamma, iInter, iIntra)
            alpha<- extractDraw(s$alpha, iInter, iIntra)
            # Calculate prediction
            for(i in 1:s$nObs){
              p  <- cost[i,2:s$nAlt]
              b  <- budget[i]
              g  <- gamma[i,2:s$nAlt]
              a0 <- alpha[i,1]
              ak <- alpha[i,2:s$nAlt]
              ph0<- phiP[i,1] 
              phk<- phiP[i,2:s$nAlt] 
              orderofV = rank(-phk)
              M = 1
              stopping = FALSE
              while(!stopping){ # nAlt
                use = orderofV<M
                #step2
                lambda_1  = b + sum(p*g*use)
                lambda_21 = ph0^(1/(1-a0))
                lambda_22 = sum(p*g*use*phk^(1/(1-a0)))
                lambda_2  = lambda_21 + lambda_22
                lambda    = (lambda_1/lambda_2)^(a0-1)
                if( M > sum(phk>lambda) || M > sum(phk>0)){
                  #step3
                  x0_1 = lambda_21*lambda_1
                  Xintra[i,1] = Xintra[i,1] + x0_1/lambda_2 # maybe X = X + ... to deal with draws
                  xk_1 = phk^(1/(1-ak))*lambda_1
                  Xintra[i,2:s$nAlt] = Xintra[i,2:s$nAlt] + use*(xk_1/lambda_2 - 1)*g  # maybe X = X + ... to deal with draws
                  stopping = TRUE
                } else M <- M + 1 # step 4
              } # end of "while"
            } # end of "for" loop over observations
          } # end of "for" loop over intra draws
          X <- X + Xintra/nIntra # CHECK IF THIS IS CORRECT!!!
        } # end of "for" loop over inter draws
        X <- X/nInter # CHECK IF THIS IS CORRECT!!!
        # Store prediction (maybe this needs to be nested in the for loops in a more clever way to address for draws)
        if(s$rawPrediction) XX[,,r] <- X else {
          Xm <- Xm + X/s$nRep
          Mm <- Mm + (X>0)/s$nRep
          Em <- Em + X*cost/s$nRep
          Xv <- Xv + apply(X  , MARGIN=2, function(v) (v-mean(v))^2)/s$nRep
          Mv <- Mv + apply(X>0, MARGIN=2, function(m) (m-mean(m))^2)/s$nRep
          Ev <- Ev + apply(X*cost, MARGIN=2, function(e) (e-mean(e))^2)/s$nRep
        }
        if(!silent){ if(r%%step2==0) cat(round(100*r/s$nRep,0), "%", sep="") else { if(r%%step1==0) cat(".") } }
      } # end of "for" loop over repetitions
    } else {
      # Without outside good
      for(r in 1:s$nRep){
        X <- 0*X
        for(iInter in 1:nInter){
          Xintra <- 0*X
          for(iIntra in 1:nIntra){
            # Extract appropiate values
            phiP <- extractDraw(s$sigma, iInter, iIntra)*epsL[[r]]
            phiP <- exp(extractDraw(s$V, iInter, iIntra) + phiP)*avail/cost
            gamma<- extractDraw(s$gamma, iInter, iIntra)
            alpha<- extractDraw(s$alpha, iInter, iIntra)
            # Calculate prediction
            for(i in 1:s$nObs){
              p  <- cost[i,]
              b  <- budget[i]
              g  <- gamma[i,]
              ak <- alpha[i,]
              phk<- phiP[i,] 
              orderofV = rank(-phk)
              M = 1
              stopping = FALSE
              while(!stopping){
                use = orderofV<M
                #step2
                lambda_1  = b + sum(p*g*use)
                lambda_2  = sum(p*g*use*phk^(1/(1-ak)))
                lambda    = (lambda_1/lambda_2)^(ak-1)
                if( M > sum(phk>lambda) || M > sum(phk>0)){
                  #step3
                  xk_1 = phk^(1/(1-ak))*lambda_1
                  Xintra[i,] = Xintra[i,] + use*(xk_1/lambda_2 - 1)*g
                  stopping = TRUE
                } else M <- M + 1 # step 4
              } # end of "while"
            } # end of "for" loop over observations
          } # end of "for" loop over intra draws
          X <- X + Xintra/nIntra
        } # end of "for" loop over inter draws
        X <- X/nInter
        # Store prediction (maybe this needs to be nested in the for loops in a more clever way to address for draws)
        if(s$rawPrediction) XX[,,r] <- X else {
          Xm <- Xm + X/s$nRep
          Mm <- Mm + (X>0)/s$nRep
          Em <- Em + X*cost/s$nRep
          Xv <- Xv + apply(X  , MARGIN=2, function(v) (v-mean(v))^2)/s$nRep
          Mv <- Mv + apply(X>0, MARGIN=2, function(m) (m-mean(m))^2)/s$nRep
          Ev <- Ev + apply(X*cost, MARGIN=2, function(e) (e-mean(e))^2)/s$nRep
        }
        if(!silent){ if(r%%step2==0) cat(round(100*r/s$nRep,0), "%", sep="") else { if(r%%step1==0) cat(".") } }
      } # end of "for" loop over repetitions
    }
    if(!silent) cat("\n")
    
    ### Prepare output
    if(s$rawPrediction) out <- XX else {
      out <- cbind(Xm, sqrt(Xv), Mm, sqrt(Mv), Em, sqrt(Ev))
      out <- apollo_insertRows(out, s$rows, NA)
      colN <- c("cont_mean", "cont_sd", "disc_mean", "disc_sd", "expe_mean", "expe_sd")
      colnames(out) <- paste( names(s$continuousChoice), rep(colN, each=s$nAlt), sep="_")
    }
    return(out)
    
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    if(!apollo_inputs$silent) apollo_print('Gradient not implemented for mdcev models')
    return(NA)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(mdcev_settings$mdcev_diagnostics(mdcev_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(mdcev_settings$mdcev_diagnostics(mdcev_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}