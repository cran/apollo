#' Calculates Cross-Nested Logit probabilities
#'
#' Calculates the probabilities of a Cross-nested Logit model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' For the model to be consistent with utility maximisation, the estimated value of the lambda parameter of all nests
#' should be between 0 and 1. Lambda parameters are inversely proportional to the correlation between the error terms of 
#' alternatives in a nest. If lambda=1,  there is no relevant correlation between the unobserved
#' utility of alternatives in that nest.
#' Alpha parameters inside \code{cnlStructure} should be between 0 and 1. Using a transformation to ensure
#' this constraint is satisfied is recommended for complex structures (e.g. logistic transformation).
#' @param cnl_settings List of inputs of the CNL model. User input is required for all settings except those with a default or marked as optional. 
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{\code{cnlNests}}: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named according to nests. The lambda at the root is fixed to 1, and therefore does not need to be defined.
#'                       \item \strong{\code{cnlStructure}}: Numeric matrix. One row per nest and one column per alternative. Each element of the matrix is the alpha parameter of that (nest, alternative) pair.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                       \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                     }
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
#'           \item \strong{\code{"gradient"}}: Not implemented.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{cnl_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}.
#'           \item \strong{\code{"report"}}: List with tree structure and choice overview.
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @export
apollo_cnl <- function(cnl_settings, functionality){
  ### Set or extract componentName
  modelType   = "CNL"
  if(is.null(cnl_settings[["componentName"]])){
    cnl_settings[["componentName"]] = ifelse(!is.null(cnl_settings[['componentName2']]),
                                             cnl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && cnl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 cnl_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, cnl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", cnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used
  if(!is.null(cnl_settings[["utilities"]])) names(cnl_settings)[which(names(cnl_settings)=="utilities")]="V"
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(cnl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load cnl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(cnl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded cnl_settings, restore the one received as argument
    if(is.null(tmp$V)           ) tmp$V            <- cnl_settings$V
    if(is.null(tmp$cnlNests)    ) tmp$cnlNests     <- cnl_settings$cnlNests
    if(is.null(tmp$cnlStructure)) tmp$cnlStructure <- cnl_settings$cnlStructure
    cnl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    cnl_settings <- apollo_preprocess(inputs = cnl_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print("No C++ optimisation available for CNL")
    cnl_settings$probs_CNL <- function(cnl_settings, all=FALSE){
      # Fix choiceVar if "raw" and choiceVar==NA
      cnl_settings$choiceNA = FALSE
      if(all(is.na(cnl_settings$choiceVar))){
        cnl_settings$choiceVar = cnl_settings$alternatives[1]
        cnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      cnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                               cnl_settings$V, cnl_settings$avail, SIMPLIFY=FALSE)
      # Extract chosen V or maximum V
      if(!all) VSubs <- Reduce('+', mapply("*", cnl_settings$Y, cnl_settings$V, SIMPLIFY=FALSE)) else VSubs <- do.call(pmax, cnl_settings$V)
      cnl_settings$V <- lapply(cnl_settings$V, "-", VSubs)
      rm(VSubs)
      # consider availabilities once before exponentiating (avoids issues if unavailable alternatives have attributes at 999)
      cnl_settings$V <- mapply('*', cnl_settings$V, cnl_settings$avail, SIMPLIFY = FALSE)
      # exponentiate utilities
      cnl_settings$V = lapply(X=cnl_settings$V, FUN=exp)
      # consider availabilities (it assumes eV and avail are in the same order)
      cnl_settings$V <- mapply('*', cnl_settings$V, cnl_settings$avail, SIMPLIFY = FALSE)
      # work out denominator for within nest probs
      denom_within = list()
      nests        = nrow(cnl_settings$cnlStructure)
      for(t in 1:nests){
        denom_within[[t]]=0
        for(j in 1:cnl_settings$nAlt){
          denom_within[[t]] = denom_within[[t]] + 
            (cnl_settings$cnlStructure[t,j]*cnl_settings$V[[cnl_settings$altnames[j]]])^(1/cnl_settings$cnlNests[[t]])
        }
      }
      # work out within-nest probs
      Pwithin = list()
      for(j in 1:cnl_settings$nAlt){
        Pwithin[[j]] = list()
        for(t in 1:nests){
          # includes a failsafe in denominator for empty nests, numerator will ensure ratio is 0 anyway
          Pwithin[[j]][[t]] = (cnl_settings$cnlStructure[t,j]*cnl_settings$V[[cnl_settings$altnames[j]]])^(1/cnl_settings$cnlNests[[t]])/(denom_within[[t]]+(denom_within[[t]]==0)) 
        }
      }
      # work out nest probs
      Pnest = list()
      denom_nest = 0
      for(t in 1:nests) denom_nest = denom_nest + denom_within[[t]]^cnl_settings$cnlNests[[t]]
      for(t in 1:nests) Pnest[[t]] = (denom_within[[t]]^cnl_settings$cnlNests[[t]])/denom_nest
      # work individual probs
      Palts = list()
      for(j in 1:cnl_settings$nAlt){
        Palts[[j]]=0
        for(t in 1:nests) Palts[[j]] = Palts[[j]] + Pnest[[t]]*Pwithin[[j]][[t]]
      }
      names(Palts) <- names(cnl_settings$V)
      if(!(all && cnl_settings$choiceNA)) Palts[["chosen"]] <- Reduce('+', mapply('*', cnl_settings$Y, Palts, SIMPLIFY=FALSE))
      if(!all) Palts <- Palts[["chosen"]]
      return(Palts)
    }
    
    cnl_settings$cnl_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
      
        # turn scalar availabilities into vectors
        for(i in 1:length(inputs$avail)) if(length(inputs$avail[[i]])==1) inputs$avail[[i]] <- rep(inputs$avail[[i]], inputs$nObs)
        
        # Construct summary table of availabilities and market share
        choicematrix = matrix(0, nrow=4, ncol=length(inputs$altnames), 
                              dimnames=list(c("Times available", "Times chosen", "Percentage chosen overall",
                                              "Percentage chosen when available"), inputs$altnames))
        choicematrix[1,] = unlist(lapply(inputs$avail, sum))
        for(j in 1:length(inputs$altnames)) choicematrix[2,j] = sum(inputs$choiceVar==inputs$altcodes[j]) # number of times each alt is chosen
        choicematrix[3,] = choicematrix[2,]/inputs$nObs*100 # market share
        choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
        choicematrix[4,!is.finite(choicematrix[4,])] <- 0
        
        if(!apollo_inputs$silent & data){
          if(any(choicematrix[4,]==0)) apollo_print("Some alternatives are never chosen in your data!", type="w")
          if(any(choicematrix[4,]>=100)) apollo_print("Some alternatives are always chosen when available!", type="w")
          #if(inputs$avail_set) apollo_print("Availability not provided (or some elements are NA). Full availability assumed.", type="i")
          apollo_print("\n")
          apollo_print(paste0('Overview of choices for ', toupper(inputs$modelType), ' model component ', 
                              ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
          print(round(choicematrix,2))
        }
      
        if(param){
          if(!apollo_inputs$silent & data) apollo_print('\n')
          if( !all(sapply(inputs$cnlNests, length)==1) ){
            inputs$cnlNests <- lapply(inputs$cnlNests, mean)
            if(!apollo_inputs$silent) apollo_print(paste0('In the following table, lambda values ',
                                                          '(nest scales) where averaged across ',
                                                          'participants and draws for model component "',
                                                          inputs$componentName, '".')) 
          }
          out_tree = cbind(inputs$cnlStructure, unlist(inputs$cnlNests))
          out_tree = apply(out_tree, MARGIN=2, function(x){
            if(all(x %in% 0:1)) round(x,0) else round(x,4)
            return(x)
          } )
          ###change 4 Oct
          #rownames(out_tree)=inputs$nestnames
          rownames(out_tree)=paste0(inputs$nestnames," nest")
          #colnames(out_tree)=c(inputs$altnames,"lambda")
          colnames(out_tree)=c(paste0(inputs$altnames," (alpha)"),"lambda")
          if(!apollo_inputs$silent){
            apollo_print(paste0('Structure for ', toupper(inputs$modelType), ' model component ', 
                                ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
            apollo_print(out_tree)
            for(a in 1:(ncol(out_tree)-1)){
              #if(sum(out_tree[,a])!=1) apollo_print(paste0("Allocation parameters for alternative \'",inputs$altnames[a],"\' do not sum to 1. You may want to impose a constraint using a logistic transform."))
              if(abs(sum(out_tree[,a])-1)>10^-4) apollo_print(paste0("Allocation parameters for alternative \'",inputs$altnames[a],"\' do not sum to 1. You may want to impose a constraint using a logistic transform."))
              if(any(out_tree[,a]<0)) apollo_print(paste0("Some allocation parameters for alternative \'",inputs$altnames[a],"\' are negative. You may want to impose a constraint using a logistic transform."))
              if(any(out_tree[,a]>1)) apollo_print(paste0("Some allocation parameters for alternative \'",inputs$altnames[a],"\' are larger than 1. You may want to impose a constraint using a logistic transform."))
            }
            if(any(out_tree[,ncol(out_tree)]<0)) apollo_print("Some lambda parameters are negative. You may want to impose a constraint or reconsider your model structure.")
            if(any(out_tree[,ncol(out_tree)]>1)) apollo_print("Some lambda parameters are larger than 1. You may want to impose a constraint or reconsider your model structure.")
          }
        }
        return(invisible(TRUE))
    }
    
    
    # Store model type
    cnl_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && functionality %in% c("preprocess", "gradient")
    test <- test && all(sapply(cnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    cnl_settings$gradient <- FALSE
    if(test){
      cnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, cnl_settings$V)
      cnl_settings$gradient <- !is.null(cnl_settings$dV)
    }; rm(test)
    
    # Return cnl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      cnl_settings$V            <- NULL
      cnl_settings$cnlNests     <- NULL
      cnl_settings$cnlStructure <- NULL
      return(cnl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(cnl_settings$V, is.function))){
    cnl_settings$V = lapply(cnl_settings$V, function(f) if(is.function(f)) f() else f )
  }
  if(any(sapply(cnl_settings$cnlNests, is.function))){
    cnl_settings$cnlNests = lapply(cnl_settings$cnlNests, function(f) if(is.function(f)) f() else f )
  }
  if(is.function(cnl_settings$cnlStructure)) cnl_settings$cnlStructure <- cnl_settings$cnlStructure()
  cnl_settings$V <- lapply(cnl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V if neccesary
  cnl_settings$V        <- cnl_settings$V[cnl_settings$altnames]
  #cnl_settings$cnlNests <- cnl_settings$cnlNests[cnl_settings$altnames]
  if(!all(cnl_settings$rows)){
    cnl_settings$V        <- lapply(cnl_settings$V, apollo_keepRows, r=cnl_settings$rows)
    #cnl_settings$cnlNests <- lapply(cnl_settings$cnlNests, apollo_keepRows, r=cnl_settings$rows)
  } 
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(cnl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) cnl_settings$cnl_diagnostics(cnl_settings, apollo_inputs)

    testL <- cnl_settings$probs_CNL(cnl_settings)
    if(any(!cnl_settings$rows)) testL <- apollo_insertRows(testL, cnl_settings$rows, 1) # insert excluded rows with value 1
    if(all(testL==0)) stop('CALCULATION ISSUE - All observations have zero probability at starting value for model component "', cnl_settings$componentName,'"')
        if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', cnl_settings$componentName,'"'))
    return(invisible(testL))
  }

  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:cnl_settings$nAlt) if(length(cnl_settings$avail[[i]])==1) cnl_settings$avail[[i]] <- rep(cnl_settings$avail[[i]], cnl_settings$nObs)
    # number of available alts in each observation
    nAvAlt <- rowSums(matrix(unlist(cnl_settings$avail), ncol=cnl_settings$nAlt))
    P = 1/nAvAlt # likelihood at zero
    if(any(!cnl_settings$rows)) P <- apollo_insertRows(P, cnl_settings$rows, 1)
    return(P)
  }

  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality=="shares_LL"){
    for(i in 1:length(cnl_settings$avail)) if(length(cnl_settings$avail[[i]])==1) cnl_settings$avail[[i]] <- rep(cnl_settings$avail[[i]], cnl_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, cnl_settings$avail)) # number of available alts in each observation
    Y = do.call(cbind,cnl_settings$Y)
    if(var(nAvAlt)==0){
      Yshares = colSums(Y)/nrow(Y)
      P = as.vector(Y%*%Yshares)
    } else {
      ## Estimate model with constants only
      mnl_ll = function(b, A, Y) as.vector(Y%*%c(b,0) - log(rowSums( A%*%exp(c(b,0)) )))
      A = do.call(cbind, cnl_settings$avail)
      b = maxLik::maxLik(mnl_ll, start=rep(0, cnl_settings$nAlt - 1), 
                         method='BFGS', finalHessian=FALSE, A=A, Y=Y)$estimate
      P = exp(mnl_ll(b, A, Y))
    }
    if(any(!cnl_settings$rows)) P <- apollo_insertRows(P, cnl_settings$rows, 1)
    return(P)
  }
  
  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #

  if(functionality %in% c("estimate","conditionals", "output", "components")){
    P <- cnl_settings$probs_CNL(cnl_settings, all=FALSE)
    if(any(!cnl_settings$rows)) P <- apollo_insertRows(P, cnl_settings$rows, 1) # insert excluded rows with value 1
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- cnl_settings$probs_CNL(cnl_settings, all=TRUE)
    if(any(!cnl_settings$rows)) P <- lapply(P, apollo_insertRows, r=cnl_settings$rows, val=NA) # insert excluded rows with value NA
    return(P)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- utils::capture.output(cnl_settings$cnl_diagnostics(cnl_settings, apollo_inputs, param=FALSE))
    P$param <- utils::capture.output(cnl_settings$cnl_diagnostics(cnl_settings, apollo_inputs, data =FALSE))
    return(P)
  }
}
