#' Calculates Fractional Nested Logit probabilities
#'
#' Calculates the probabilities of a Fractional Nested Logit (FNL) model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' In this implementation of the Nested Logit model, each nest must have a lambda parameter associated to it.
#' For the model to be consistent with utility maximisation, the estimated value of the Lambda parameter of all nests
#' should be between 0 and 1. Lambda parameters are inversely proportional to the correlation between the error terms of 
#' alternatives in a nest. If lambda=1, then there is no relevant correlation between the unobserved
#' utility of alternatives in that nest.
#' The tree must contain an upper nest called \code{"root"}. The lambda parameter of the root is automatically
#' set to 1 if not specified in \code{nlNests}, but can be changed by the user if desired (though not advised).
#' @param fnl_settings List of inputs of the FNL model. It should contain the following.
#'                    \itemize{
#'                       \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'utilities'.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{choiceShares}}: Named list of numeric vectors. Share allocated to each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                       \item \strong{\code{nlNests}}: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named with the nest name. The lambda at the root is automatically fixed to 1 if not provided by the user.
#'                       \item \strong{\code{nlStructure}}: Named list of character vectors. As many elements as nests, it must include the "root". Each element contains the names of the nests or alternatives that belong to it. Element names must match those in \code{nlNests}.
#'                       \item \strong{\code{utilities}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                    }
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
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{fnl_settings}.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: List with tree structure and choice overview.
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @importFrom utils capture.output
#' @export
apollo_fnl <- function(fnl_settings, functionality){
  ### Set or extract componentName
  modelType   = "FNL"
  if(is.null(fnl_settings[["componentName"]])){
    fnl_settings[["componentName"]] = ifelse(!is.null(fnl_settings[['componentName2']]),
                                            fnl_settings[['componentName2']], modelType)
    test <- functionality=="validate" && fnl_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 fnl_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, fnl_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", fnl_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  #### replace utilities by V if used
  if(!is.null(fnl_settings[["utilities"]])) names(fnl_settings)[which(names(fnl_settings)=="utilities")]="V"
  
  ### turn character vector of alternatives into named vector (to make compatible with MNL code)
  if(is.character(fnl_settings[["alternatives"]])) fnl_settings[["alternatives"]]=setNames(1:length(fnl_settings[["alternatives"]]), fnl_settings[["alternatives"]])
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(fnl_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load fnl_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(fnl_settings$componentName, "_settings")]]
    # If there is no V inside the loaded fnl_settings, restore the one received as argument
    if(is.null(tmp$V)          ) tmp$V           <- fnl_settings$V
    if(is.null(tmp$nlNests)    ) tmp$nlNests     <- fnl_settings$nlNests
    if(is.null(tmp$nlStructure)) tmp$nlStructure <- fnl_settings$nlStructure
    fnl_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    fnl_settings <- apollo_preprocess(inputs = fnl_settings, modelType, 
                                     functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print("No C++ optimisation available for FNL")
    fnl_settings$probs_FNL <- function(fnl_settings, all=FALSE){
      # Fix choiceVar if "raw" and choiceVar==NA
      fnl_settings$choiceNA = FALSE
      if(all(is.na(unlist(fnl_settings$choiceShares)))){
        for(s in 1:length(fnl_settings$choiceShares)) fnl_settings$choiceShares[[s]] = rep(1/length(fnl_settings$choiceShares),length(fnl_settings$alternatives[1]))
        fnl_settings$choiceNA = TRUE
      }
      # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
      fnl_settings$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), 
                              fnl_settings$V, fnl_settings$avail, SIMPLIFY=FALSE)
      # Extract chosen V or maximum V
      if(!all){
        VSubs <- mapply("*", fnl_settings$maxShare, fnl_settings$V, SIMPLIFY=FALSE)
        VSubs <- Reduce('+', VSubs)
      }else{
        VSubs <- do.call(pmax, fnl_settings$V)
      } 
      fnl_settings$V <- lapply(fnl_settings$V, "-", VSubs)
      rm(VSubs)
      # Not sure what the two following lines are supposed to used for
      #combined_elements="root"
      #for(j in 1:length(nlStructure)) combined_elements=c(combined_elements,nlStructure[[j]])
      # loop over nests to create new utility elements and new availability terms
      for(k in length(fnl_settings$nlStructure):1){
        nestK <- names(fnl_settings$nlStructure)[k]
        fnl_settings$V[[nestK]] = 0
        # calculate availability of nest
        fnl_settings$avail[[nestK]] = 1*( Reduce('+', fnl_settings$avail[ fnl_settings$nlStructure[[k]] ]) > 0 )
        for(j in 1:length(fnl_settings$nlStructure[[k]])){
          nodeJ <- fnl_settings$nlStructure[[k]][j]
          fnl_settings$V[[nestK]] = fnl_settings$V[[nestK]] + 
            fnl_settings$avail[[nodeJ]]*exp( fnl_settings$V[[nodeJ]]/fnl_settings$nlNests[[nestK]] )
        }
        fnl_settings$V[[nestK]] = fnl_settings$nlNests[[nestK]]*log(fnl_settings$V[[nestK]])
      }
      # calculate log(probabilities)
      logPalts=list()
      for(j in 1:length(fnl_settings$altnames)){
        logPalts[[j]]=0
        ancestorsJ <- fnl_settings$ancestors[[fnl_settings$altnames[[j]]]]
        for(k in 1:(length(ancestorsJ)-1)){ # loop to level just below root
          current_V = fnl_settings$V[[ ancestorsJ[k] ]]
          next_V    = fnl_settings$V[[ ancestorsJ[k+1] ]]
          logPalts[[j]] = logPalts[[j]] + (current_V-next_V)/fnl_settings$nlNests[[ ancestorsJ[k+1] ]]
        }
      }
      Palts = lapply(X=logPalts, FUN=exp)
      names(Palts)=names(fnl_settings$V)[1:length(fnl_settings$altnames)]
      # consider availabilities (it assumes Palts and avail are in the same order)
      Palts <- mapply('*', Palts, fnl_settings$avail[1:length(fnl_settings$altnames)], SIMPLIFY = FALSE)
      Palts <- lapply(Palts, function(x) {
        x[is.na(x)] <- 0
        return(x)}) # replace all NaN by 0
      # Prepare output
      if(!(all && fnl_settings$choiceNA)){
        Palts[["chosen"]] <- Reduce('*', mapply('^', Palts, fnl_settings$choiceShares, SIMPLIFY=FALSE))
      } 
      if(!all) Palts <- Palts[["chosen"]]
      return(Palts)
    }
    
    fnl_settings$fnl_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
      
      # Table describing dependent variable
      choicematrix <- matrix(0, nrow=4, ncol=inputs$nAlt, 
                             dimnames=list(c("Times available",
                                             "Observations with non-zero share",
                                             "Average share overall",
                                             "Average share when available"),
                                           inputs$altnames))
      for(a in 1:inputs$nAlt){
        choicematrix[1,a] <- ifelse(length(inputs$avail[[a]])==1 && inputs$avail[[a]]==1, 
                                    inputs$nObs, sum(inputs$avail[[a]]) )
        choicematrix[2,a] <- sum(inputs$choiceShares[[a]]>0)
        choicematrix[3,a] <- ifelse(choicematrix[1,a]>0, mean(inputs$choiceShares[[a]]), 0)
        choicematrix[4,a] <- ifelse(choicematrix[1,a]>0, sum(inputs$choiceShares[[a]])/choicematrix[1,a], 0)
      }
      # Print table
      if(!apollo_inputs$silent & data){
        apollo_print("\n")
        apollo_print(paste0('Overview of choices for ', toupper(inputs$modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        print(round(choicematrix,2))
        
        # Print warnings
        for(a in 1:inputs$nAlt){
          if(choicematrix[2,a]==0) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'), type="w")
          if(choicematrix[4,a]==1) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always given the full choice share when available in model component "', inputs$componentName, '".'), type="w")
        }
        #if(inputs$avail_set==TRUE & !apollo_inputs$silent) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'), type="w")
      }
      

        if(param){
          if(!apollo_inputs$silent & data) apollo_print('\n') # 
          if(!apollo_inputs$silent){
            # WARNING for automatic setting of root nesting parameter
            if(inputs$root_set) apollo_print("Notice: Root lambda parameter set to 1.")
            # Identifying nest's parents
            nestAbove <- unique(lapply(inputs$ancestors, '[', -1))
            nestAbove <- setNames(sapply(nestAbove, function(x) if(length(x)==1) return('Inf') else x[2]) ,
                                  sapply(nestAbove, '[', 1))
            # Printing graphical representation of the tree, using recursive function
            apollo_print(paste0('Nesting structure for ', toupper(inputs$modelType), ' model component ', 
                                ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
            print_tree_level = function(nlStructure, component, preceding_nest_layer, space){
              if(preceding_nest_layer!=0) space=c(space,"  |")
              for(j in 1:length(nlStructure[[component]])){
                space <- gsub("[']", " ", space)
                if(j==length(nlStructure[[component]])) space[length(space)] <- gsub("[|]", "'", space[length(space)])
                if(nlStructure[[component]][j] %in% inputs$altnames){
                  depth <- length(space)
                  cat("\n",space,rep("-",3*(maxDepth-depth)),"-Alternative: ",nlStructure[[component]][j], sep="")
                } else {
                  l  <- inputs$nlNests[[nlStructure[[component]][j]]]
                  #n0 <- nestAbove[nlStructure[[component]][j]]
                  #if(n0=='Inf') l0 <- 1 else l0 <- inputs$nlNests[[n0]]
                  if(length(l)>1){
                    cat("\n",space,"-Nest: ", nlStructure[[component]][j], " (distributed, mean: ",mean(l),")", sep="")  
                  }else{
                    cat("\n",space,"-Nest: ", nlStructure[[component]][j], " (",round(l,4), ")", sep="")
                  }
                  #if(any(l<0 | l0<l)) cat(' WARNING: nest param. should be between 0 and ', round(l0,4), '.', sep='')
                  print_tree_level(nlStructure, nlStructure[[component]][j], preceding_nest_layer+1, space)
                }
              }
            } # end of print_tree_level function
            maxDepth <- max(sapply(inputs$ancestors, length))-1
            cat("Nest: ",names(inputs$nlStructure)[[1]]," (",round(inputs$nlNests[[names(inputs$nlStructure)[[1]]]],4),")", sep="")
            print_tree_level(inputs$nlStructure, "root", preceding_nest_layer=0, space="|")
            apollo_print('\n')
            # Print warning if nesting parameters do not make sense
            for(i in names(inputs$nlNests)){
              l  <- inputs$nlNests[[i]]
              if(i=='root') l0 <- 1 else l0 <- inputs$nlNests[[ nestAbove[i] ]]
              #if(any(l<0 | l0<l)){
              if(length(l)==1 && any(l<0 | l0<l)){
                txt <- paste0('The nesting parameter for nest "', i, '" should be between 0 and ', round(l0,4))
                if(i!='root') txt <- paste0(txt, ' (the nesting parameter for nest "', nestAbove[i], '")')
                txt <- paste0(txt, ', yet its value is ', round(l, 4), '.')
                cat('\n'); apollo_print(txt, type="w")
              }
            }
          }
        } # end of FNL special checks
      
        return(invisible(TRUE))
    }
    
    
    # Store model type
    fnl_settings$modelType <- modelType
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && functionality %in% c("preprocess", "gradient")
    test <- test && all(sapply(fnl_settings$V, is.function))
    test <- test && apollo_inputs$apollo_control$analyticGrad
    fnl_settings$gradient <- FALSE
    if(test){
      fnl_settings$dV       <- apollo_dVdB(apollo_beta, apollo_inputs, fnl_settings$V)
      fnl_settings$gradient <- !is.null(fnl_settings$dV)
    }; rm(test)
    
    # Return fnl_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      fnl_settings$V           <- NULL
      fnl_settings$nlNests     <- NULL
      fnl_settings$nlStructure <- NULL
      return(fnl_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(fnl_settings$V, is.function))){
    fnl_settings$V = lapply(fnl_settings$V, function(f) if(is.function(f)) f() else f )
  }
  if(any(sapply(fnl_settings$nlNests, is.function))){
    fnl_settings$nlNests = lapply(fnl_settings$nlNests, function(f) if(is.function(f)) f() else f )
  }
  if(is.function(fnl_settings$nlStructure)) fnl_settings$nlStructure <- fnl_settings$nlStructure()
  fnl_settings$V <- lapply(fnl_settings$V, function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  
  ### Reorder V if necessary
  fnl_settings$V        <- fnl_settings$V[fnl_settings$altnames]
  if(!all(fnl_settings$rows)) fnl_settings$V <- lapply(fnl_settings$V, apollo_keepRows, r=fnl_settings$rows)
  # No need to drop rows in avail, choiceVar nor Y, as these are
  # already filtered due to them not changing across iterations.
  
#  if(fnl_settings$root_set) fnl_settings$nlNests$root=1

  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(fnl_settings, modelType, 
                                                                   functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) fnl_settings$fnl_diagnostics(fnl_settings, apollo_inputs)
    
    testL=fnl_settings$probs_FNL(fnl_settings)
    if(any(!fnl_settings$rows)) testL <- apollo_insertRows(testL, fnl_settings$rows, 1) # insert excluded rows with value 1
    if(all(testL==0)) stop('CALCULATION ISSUE - All observations have zero probability at starting value for model component "', fnl_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', fnl_settings$componentName,'"'), type="i")
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:fnl_settings$nAlt) if(length(fnl_settings$avail[[i]])==1) fnl_settings$avail[[i]] <- rep(fnl_settings$avail[[i]], fnl_settings$nObs)
    # number of available alts in each observation
    nAvAlt <- rowSums(matrix(unlist(fnl_settings$avail), ncol=fnl_settings$nAlt))
    P = 1/nAvAlt # likelihood at zero
    if(any(!fnl_settings$rows)) P <- apollo_insertRows(P, fnl_settings$rows, 1)
    return(P)
  }
  
  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality=="shares_LL"){
    return(NA)
  }
  
  # ############################################################################ #
  #### functionality="estimate/prediction/conditionals/raw/output/components" ####
  # ############################################################################ #
  
  if(functionality %in% c("estimate","conditionals", "output", "components")){
    P <- fnl_settings$probs_FNL(fnl_settings, all=FALSE)
    if(any(!fnl_settings$rows)) P <- apollo_insertRows(P, fnl_settings$rows, 1) # insert excluded rows with value 1
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- fnl_settings$probs_FNL(fnl_settings, all=TRUE)
    if(any(!fnl_settings$rows)) P <- lapply(P, apollo_insertRows, r=fnl_settings$rows, val=NA) # insert excluded rows with value 1
    return(P)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(fnl_settings$fnl_diagnostics(fnl_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(fnl_settings$fnl_diagnostics(fnl_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
