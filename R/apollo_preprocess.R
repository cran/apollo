#' Pre-process input for multiple models
#' return
#' @param inputs List of settings
#' @param modelType Character. Type of model, e.g. "mnl", "nl", "cnl", etc.
#' @param functionality Character. Either "estimate","prediction","validate","zero_LL","conditionals","output","raw", or "preprocess". Only used for validation, it does not influence the return values.
#' @param apollo_inputs List of main inputs to the model estimation process. See \link{apollo_validateInputs}.
#' @return The returned object depends on the value of argument operation
#' @export
apollo_preprocess <- function(inputs, modelType, functionality, apollo_inputs){
  inputs$modelType <- modelType
  modelType <- tolower(modelType)
  
  #### Checks common to all models ####
  # inputs
  if(is.null(names(inputs))) stop('All elements inside the inputs lists for model components must be named, e.g. mnl_settings=list(alternatives=c(...), avail=...).')
  if(is.null(inputs[["componentName"]])) stop('The settings of at least one model component is missing the mandatory "componentName" object.')
  # functionality
  test <- functionality %in% c("estimate","prediction","validate","zero_LL","conditionals","output","raw","preprocess", "components", "gradient", "report")
  if(!test) stop("Non-permissable setting for \"functionality\" for model component \"",inputs$componentName,"\"")
  
  #### MNL, NL, CNL, DFT ####
  if(modelType %in% c("mnl","nl","cnl","dft")){
    # Check for mandatory inputs
    mandatory <- c("alternatives", "choiceVar")
    if(modelType!="dft") mandatory <- c(mandatory, "V")
    if(modelType=="cnl") mandatory <- c(mandatory, "cnlNests", "cnlStructure")
    if(modelType=="nl") mandatory <- c(mandatory, "nlNests", "nlStructure")
    if(modelType=="dft") mandatory <- c(mandatory, "attrValues", "altStart", "attrWeights", "attrScalings", "procPars")
    for(i in mandatory) if(!(i %in% names(inputs))) stop('The inputs list for model component "', inputs$componentName, 
                                                         '" needs to include an object called "', i,'"!')
    # Check for optional inputs (avail and rows)
    if(is.null(inputs[["rows"]])) inputs[["rows"]]="all"
    if(is.null(inputs[['avail']])){
      inputs[['avail']]=1
      if(!apollo_inputs$silent && functionality=='validate') apollo_print('Setting "avail" is missing, so full availability is assumed.')
    }
    
    ### Store useful values
    inputs$altnames = names(inputs$alternatives)
    inputs$altcodes = inputs$alternatives
    inputs$nAlt     = length(inputs$alternatives)
    inputs$nObs <- tryCatch(if(!is.null(apollo_inputs$database)) nrow(apollo_inputs$database) else stop('x'),
                            error=function(e){
                              lenV <- sapply(inputs$V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                              lenA <- sapply(inputs$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                              lenC <- length(inputs$choiceVar)
                              return(max(lenV, lenA, lenC))
                            })
    if(modelType=='cnl') inputs$nestnames <- names(inputs$cnlNests)
    if(modelType=='nl' ) inputs$nestnames <- names(inputs$nlStructure)
    if(modelType=='dft'){
      mandatory2 <- c("error_sd", "timesteps", "phi1", "phi2")
      for(i in mandatory2) if(!(i %in% names(inputs$procPars))) stop('The inputs list for model component "', inputs$componentName, 
                                                                     '" needs to include an object called "procPars" with an element called "', i,'"!')
      inputs$nAttrs <- max(length(inputs$attrWeights), length(inputs$attrScalings))
      inputs$warn1 <- FALSE
      inputs$warn2 <- FALSE
      inputs$warn3 <- FALSE
      inputs$Dims = 1
      if(!is.null(apollo_inputs$apollo_control$mixing) && apollo_inputs$apollo_control$mixing){
        if(apollo_inputs$apollo_draws$interNDraws!=0) inputs$Dims = 2
        if(apollo_inputs$apollo_draws$intraNDraws!=0) inputs$Dims = 3
      }
    }
    
    ### Format checks
    # alternatives
    test <- is.vector(inputs$alternatives) & !is.null(names(inputs$alternatives))
    if(!test) stop("The \"alternatives\" argument for model component \"",inputs$componentName,"\" needs to be a named vector")
    # avail
    test <- is.list(inputs$avail) || (length(inputs$avail)==1 && inputs$avail==1)
    if(!test) stop("The \"avail\" argument for model component \"",inputs$componentName,"\" needs to be a list or set to 1")
    if(is.list(inputs$avail)){
      lenA <- sapply(inputs$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
      test <- all(lenA==inputs$nObs | lenA==1)
      if(!test) stop("All entries in \"avail\" for model component \"",inputs$componentName,"\" need to be a scalar or a vector with one entry per observation in the \"database\"")
    }
    # choiceVar
    test <- is.vector(inputs$choiceVar) && (length(inputs$choiceVar)==inputs$nObs || length(inputs$choiceVar)==1)
    if(!test) stop("The \"choiceVar\" argument for model component \"",inputs$componentName,"\" needs to be a scalar or a vector with one entry per observation in the \"database\"")
    if(modelType!="dft"){
      # V
      if(!is.list(inputs$V)) stop("The \"V\" argument for model component \"",inputs$componentName,"\" needs to be a list")
      lenV <- sapply(inputs$V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
      test <- all(lenV==inputs$nObs | lenV==1)
      if(!test) stop("Each element of \"V\" for model component \"",inputs$componentName,"\" needs to be a scalar or a vector/matrix/cube with one row per observation in the \"database\"")  
    }
    # rows
    test <- is.vector(inputs$rows)
    test <- test && ( (is.logical(inputs$rows) && length(inputs$rows)==inputs$nObs) || (length(inputs$rows)==1 && inputs$rows=="all") )
    if(!test) stop("The \"rows\" argument for model component \"",inputs$componentName,"\" needs to be \"all\" or a vector of logical statements with one entry per observation in the \"database\"")
    
    if(modelType=="cnl"){
      # cnlNests
      test <- is.list(inputs$cnlNests) && !is.null(names(inputs$cnlNests)) && 
        all(sapply(inputs$cnlNests, is.vector)) && all(sapply(inputs$cnlNests, is.numeric))
      if(!test) stop('Argument "cnlNests" for model component "', inputs$componentName,'" must be a named list of numeric vectors describing the lambda parameter of each nest.')
      len  <- sapply(inputs$cnlNests, function(x) ifelse(is.array(x), dim(x)[1], length(x)) )
      test <- all(len==inputs$nObs | len==1)
      if(!test) stop('Each element of "cnlNests" for model component "', inputs$componentName,'" must be a scalar or a vector with one entry per observation in the "database"')
      # cnlStructure
      test <- is.matrix(inputs$cnlStructure) && is.numeric(inputs$cnlStructure)
      if(!test) stop('Argument "cnlStructure" for model component "', inputs$componentName,'" must be a numeric matrix.')
    }
    
    if(modelType=='nl'){
      # nlStructure
      test <- is.list(inputs$nlStructure) && !is.null(names(inputs$nlStructure)) 
      test <- test && all(sapply(inputs$nlStructure, is.vector)) && all(sapply(inputs$nlStructure, is.character))
      if(!test) stop('Argument "nlStructure" for model component "', inputs$componentName,'" must be a named list of character vectors describing the contents of each nest.')
      # Set lambda_root to 1 if necessary
      inputs$root_set  <- FALSE
      if(!("root" %in% names(inputs$nlNests))){
        inputs$root_set <- TRUE
        inputs$nlNests["root"] <- 1
      }
      # Order tree structure
      nlStructure_ordered=list()
      element_list="root"
      j=1
      while(j>0){
        temp = rep(TRUE,length(element_list))
        for(k in 1:length(element_list)) if(element_list[k] %in% inputs$altnames) temp[k]=FALSE
        element_list = element_list[temp]
        if(length(element_list)>0){
          nlStructure_ordered[[element_list[1]]] = inputs$nlStructure[[element_list[1]]]
          element_list = c(element_list, inputs$nlStructure[[element_list[1]]])
          element_list = element_list[-1]
        }
        j = length(element_list)
      }
      inputs$nlStructure = nlStructure_ordered
      # Calculate ancestors
      ancestors=list()
      for(j in 1:length(inputs$altnames)){
        altJ <- inputs$altnames[[j]]
        ancestors[[altJ]] = altJ
        current = altJ
        for(k in length(inputs$nlStructure):1){
          if(current %in% inputs$nlStructure[[k]]){
            ancestors[[inputs$altnames[[j]]]] = c(ancestors[[altJ]], names(inputs$nlStructure)[k])
            current = names(inputs$nlStructure)[k]
          }
        }
      }
      inputs$ancestors <- ancestors
    }
    
    ### Expand availabilities if necessary
    inputs$avail_set <- FALSE
    if(length(inputs$avail)==1 && inputs$avail==1){
      inputs$avail <- as.list(setNames(rep(1,inputs$nAlt), inputs$altnames))
      inputs$avail_set <- TRUE
    }
    
    ### Check that avail and V are available for all alternatives
    if(modelType!="dft" && !all(inputs$altnames %in% names(inputs$V))) stop("The names of the alternatives for model component \"", inputs$componentName,"\" do not match those in \"V\".")
    if(!all(inputs$altnames %in% names(inputs$avail))) stop("The names of the alternatives for model component \"", inputs$componentName,"\" do not match those in \"avail\".")
    
    ### Reorder availabilities
    inputs$avail <- inputs$avail[inputs$altnames]
    
    ### Expand rows if necessary, and update nObs
    if(length(inputs$rows)==1 && inputs$rows=="all") inputs$rows <- rep(TRUE, inputs$nObs)
    inputs$nObs <- sum(inputs$rows)
    # Filter rows, except for V
    if(any(!inputs$rows)){
      inputs$avail <- lapply(inputs$avail, 
                             function(av) if(length(av)==1) return(av) else return(av[inputs$rows]))
      inputs$choiceVar <- apollo_keepRows(inputs$choiceVar, inputs$rows)
      if(modelType=="cnl") inputs$cnlNests <- lapply(inputs$cnlNests, 
                                                     function(x) if(length(x)==1) return(x) else return(x[inputs$rows]))
    }
    
    if(modelType=="dft"){
      #### check that either attrWeights or attrScalings is supplied, but not both
      s1 = sum(lengths(inputs[["attrWeights"]]))
      if(s1>1) inputs$attrnames = names(inputs[["attrWeights"]]) else inputs$attrnames = names(inputs[['attrScalings']])
      
      # Check that the elements of attrValues match altnames
      test <- all(names(inputs$attrValues) %in% inputs$altnames)
      if(!test) stop('The "attrValues" attribute names for model component "', inputs$componentName,
                     '" do not match those given in "alternatives"!') 
      
      # Check for additional warnings
      # Give warning message if any of the elements in attrValues are not in attrnames
      for (i in 1:inputs$nAlt) if(!all(names(inputs$attrValues[[i]]) %in% inputs$attrnames)) inputs$warn1 <- TRUE 
      # Give warning message if any of the elements in altStarts are not in altnames
      if(!all(names(inputs$altStart) %in% inputs$altnames)) inputs$warn2 <- TRUE 
      # Check altStart is a list:
      if(!is.list(inputs$altStart)) {
        inputs$altStart=list()
        inputs$warn3 <- TRUE
      }
      
      # Add zeros to attrValues for attributes not supplied.
      for(i in 1:inputs$nAlt) for(j in 1:inputs$nAttrs){
        test <- is.null(inputs$attrValues[[inputs$altnames[i]]][[inputs$attrnames[j]]])
        if(test) inputs$attrValues[[inputs$altnames[i]]][[inputs$attrnames[j]]]=0
      }
      
      ## check that each alternative attribute is of length nObs and update if necessary (will extend zeros, for examples)
      for (i in 1:inputs$nAttrs) for (j in 1:inputs$nAlt) if(length(inputs$attrValues[[j]][[i]])==1){
        inputs$attrValues[[j]][[i]] = rep(inputs$attrValues[[j]][[i]],inputs$nObs)
      }
    }
    
    ### Create Y
    inputs$Y <- lapply(as.list(inputs$alternatives), function(i) inputs$choiceVar==i)
    
    # Record availability of chosen alternative
    inputs$chosenAvail <- Reduce('+', mapply('*', inputs$Y, inputs$avail, SIMPLIFY=FALSE))
  }
  
  #### Exploded logit ####
  if(modelType=="el"){
    # Check for mandatory inputs
    mandatory <- c("alternatives", "choiceVars", "V")
    for(i in mandatory) if(!(i %in% names(inputs))) stop('The inputs list for model component "', inputs$componentName, 
                                                         '" needs to include an object called "', i,'"!')
    # Check for optional inputs
    if(is.null(inputs[["rows"]]  )) inputs[["rows"]]="all"
    if(is.null(inputs[["scales"]])){
      inputs[["scales"]] <- as.list( rep(1, length(inputs$choiceVars)) ) # Changed length to number of choices (to allow for incomplete explosion) 7/05/2020
      inputs$fixedScales <- TRUE
    } else inputs$fixedScales <- FALSE
    if(is.null(inputs[['avail']])){
      inputs[['avail']]=1
      if(!apollo_inputs$silent && functionality=='validate') apollo_print('Setting "avail" is missing, so full availability is assumed.')
    }
    
    ### Store useful values
    inputs$altnames = names(inputs$alternatives)
    inputs$altcodes = inputs$alternatives
    inputs$nAlt     = length(inputs$alternatives)
    inputs$nObs <- tryCatch(nrow(apollo_inputs$database),
                            error=function(e){
                              lenV <- sapply(inputs$V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                              lenA <- sapply(inputs$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                              lenC <- length(inputs$choiceVar)
                              return(max(lenV, lenA, lenC))
                            })
    inputs$stages = length(inputs$choiceVars)
    
    # Expand availabilities if necessary
    inputs$avail_set <- FALSE
    if(length(inputs$avail)==1 && inputs$avail==1){
      inputs$avail <- as.list(setNames(rep(1,inputs$nAlt), inputs$altnames))
      inputs$avail_set <- TRUE
    }
    
    # alternatives
    test <- is.vector(inputs$alternatives) & !is.null(names(inputs$alternatives))
    if(!test) stop("The \"alternatives\" argument for model component \"",inputs$componentName,"\" needs to be a named vector")
    if(-1 %in% inputs$altcodes) stop("Negative one (-1) is not valid code for an alternative!") # Because it is used to represent a lack of a choice
    
    # avail
    if(!is.list(inputs$avail) && !is.null(names(inputs$avail)) ) stop("The \"avail\" argument for model component \"",inputs$componentName,"\" needs to be a named list or set to 1")
    lenA <- sapply(inputs$avail, function(v) if(is.array(v)) dim(v)[1] else length(v) )
    test <- all(lenA==inputs$nObs | lenA==1)
    if(!test) stop("All entries in \"avail\" for model component \"",inputs$componentName,"\" need to be a scalar or a vector with one entry per observation in the \"database\"")
    if(!all(inputs$altnames %in% names(inputs$avail))) stop('The names of the alternatives for model component "',
                                                            inputs$componentName,'" do not match those in "avail".')
    inputs$avail <- inputs$avail[inputs$altnames]
    
    # choiceVars
    if(!is.list(inputs$choiceVars)) stop('The "choiceVars" argument for model component "',inputs$componentName,'" needs to be a list, with one vector per stage, each with one entry per observation in the "database"')
    if(length(inputs$choiceVars)<2) stop('The "choiceVars" argument for model component "',inputs$componentName,'" needs to be a list, with at least two choice variables (i.e. two choice stages).')
    test <- all(sapply(inputs$choiceVars, function(cho) is.vector(cho) && (length(cho) %in% c(1, inputs$nObs)) ))
    if(!test) stop('Each element of the "choiceVars" argument for model component "', inputs$componentName, '" should be a vector of length one or nObs.')
    inputs$choiceVars = lapply(inputs$choiceVars, function(l){l[is.na(l)] <- -1; return(l)})
    test <- all(unique(unlist(inputs$choiceVars)) %in% c(inputs$altcodes,-1))
    if(!test) stop("The data contains values in \"choiceVars\" for model component \"",inputs$componentName,"\" that are not included in \"alternatives\".")
    cho <- do.call(cbind, inputs$choiceVars) # Check nothing is chosen twice
    for(i in inputs$altcodes){
      tmp <- which( apply(cho, MARGIN=1, function(r) sum(r==inputs$altcodes[i]))>1 )
      if(length(tmp)>1) tmp <- paste0(tmp, collapse=", ")
      if(length(tmp)>0) stop('Alternative ', inputs$altcodes[i], ' chosen more than once in row(s) ', 
                             tmp, ' for model component "', inputs$componentName)
    }
    
    # V
    if(!is.list(inputs$V) || is.null(names(inputs$V))) stop('The "V" argument for model component "', inputs$componentName, '" needs to be a named list.')
    lenV <- sapply(inputs$V, function(v) if(is.array(v)) dim(v)[1] else length(v) )
    test <- all(lenV==inputs$nObs | lenV==1)
    if(!test) stop("Each element of \"V\" for model component \"",inputs$componentName,"\" needs to be a scalar or a vector/matrix/cube with one row per observation in the \"database\"")
    if(!all(inputs$altnames %in% names(inputs$V))) stop('The names of the alternatives for model component "',inputs$componentName,'" do not match those in "V".')
    
    # rows
    test <- is.vector(inputs$rows) && length(inputs$rows)==1 && inputs$rows=="all"
    test <- test || ( is.vector(inputs$rows) && is.logical(inputs$rows) && length(inputs$rows)==inputs$nObs  )
    if(!test) stop("The \"rows\" argument for model component \"",inputs$componentName,"\" needs to be \"all\" or a vector of logical statements with one entry per observation in the \"database\"")
    if(length(inputs$rows)==1 && inputs$rows=="all") inputs$rows <- rep(TRUE, inputs$nObs)
    
    # scales
    test <- is.list(inputs$scales) && length(inputs$scales)==inputs$stages
    if(!test) stop("The object \"scales\" for model component \"",inputs$componentName,"\" needs to be the same length as the number of stages!")
    
    # Create Y
    inputs$Y <- vector(mode="list", length=inputs$stages)
    for(s in 1:inputs$stages) for(a in 1:length(inputs$alternatives)) inputs$Y[[s]][[a]] <- inputs$choiceVars[[s]]==inputs$alternatives[[a]]
    for(s in 1:inputs$stages) names(inputs$Y[[s]]) <- names(inputs$alternatives)
    
    # Create new availabilities list, with one list per stage
    avail2 <- list(inputs$avail)
    for(s in 2:inputs$stages) avail2[[s]] = mapply(function(a,y) a*!y, avail2[[s-1]], inputs$Y[[s-1]], SIMPLIFY=FALSE)
    inputs$avail <- avail2
    rm(avail2)
    
    # Filter rows in choiceVars, avail and Y, and update nObs
    if(any(!inputs$rows)) for(s in 1:inputs$stages){
      inputs$avail[[s]]      <- lapply(inputs$avail[[s]], apollo_keepRows, r=inputs$rows)
      inputs$choiceVars[[s]] <- apollo_keepRows(inputs$choiceVars[[s]], apollo_keepRows)
      inputs$Y[[s]]          <- lapply(inputs$Y[[s]], apollo_keepRows, r=inputs$rows)
    }
    inputs$nObs <- sum(inputs$rows)
    
    # Record availability of chosen alternative
    inputs$chosenAvail <- list(mode="list", length=inputs$stages)
    for(s in 1:inputs$stages) inputs$chosenAvail[[s]] <- Reduce("+", mapply("*", inputs$avail[[s]], inputs$Y[[s]], SIMPLIFY=FALSE))
    
  }
  
  #### NormD ####
  if(modelType=="normd"){
    if(is.null(inputs$outcomeNormal)) stop("The normalDensity_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"outcomeNormal\"!")
    if(is.null(inputs$xNormal)      ) stop("The normalDensity_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"xNormal\"!")
    if(is.null(inputs$mu)           ) stop("The normalDensity_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"mu\"!")
    if(is.null(inputs$sigma)        ) stop("The normalDensity_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"sigma\"!")
    if(is.null(inputs$rows)         ) inputs[["rows"]] <- "all"
    
    # Expand rows if necessary
    if(length(inputs$rows)==1 && inputs$rows=="all") inputs$rows <- rep(TRUE, length(inputs$outcomeNormal))
  }
  
  #### OL, OP ####
  if(modelType %in% c("ol", "op")){
    if(is.null(inputs[["outcomeOrdered"]])) stop("The ol_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"outcomeOrdered\"!")
    if(is.null(inputs[["V"]])             ) stop("The ol_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"V\"!")
    if(is.null(inputs[["tau"]])           ) stop("The ol_settings list for model component \"",inputs$componentName,"\" needs to include an object called \"tau\"!")
    if(is.null(inputs[["coding"]])        ) inputs[["coding"]] = NULL
    if(is.null(inputs[["rows"]])          ) inputs[["rows"]]   = "all"
    inputs$nObs <- tryCatch(nrow(apollo_inputs$database),
                            error=function(e){
                              lenV <- sapply(inputs$V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                              lenC <- length(inputs$outcomeOrdered)
                              return(max(lenV, lenC))
                            })
    ### Format of tau
    test <- is.vector(inputs$tau) | is.matrix(inputs$tau) | is.list(inputs$tau)
    if(!test) stop("Thresholds for Ordered Logit for model component \"",inputs$componentName,"\" need to be a vector, matrix, or list!")
    if(is.vector(inputs$tau)){
      if(length(inputs$tau)==1) stop("If provided as scalars, need at least two thresholds for Ordered Logit for model component \"",inputs$componentName,"\"!")
      inputs$tau = as.list(inputs$tau)
    } else if(is.matrix(inputs$tau)){
      if(nrow(inputs$tau)!=inputs$nObs) stop("If provided as a matrix, need one value per observation in the data for each threshold for Ordered Logit for model component \"",inputs$componentName,"\"!")
      inputs$tau = split(inputs$tau, rep(1:ncol(inputs$tau), each=nrow(inputs$tau)))
    } else if(is.list(inputs$tau)){
      if(length(inputs$tau)==1) stop("If provided as a list, the list of thresholds need at least two elements for Ordered Logit component \"",inputs$componentName,"\"!")
      if(any(sapply(inputs$tau,is.list))) stop("If provided as a list, elements in list of thresholds for Ordered Logit for model component \"",inputs$componentName,"\" cannot be lists themselves!")
      mixing <- tryCatch(apollo_inputs$apollo_control$mixing, error = function(e) return(NA))
      #if(any(sapply(inputs$tau, is.function))){
      #  tauE <- lapply(inputs$tau, function(tt) if(is.function(tt)) tt() else tt) # might be wrong environment
      #} else tauE <- inputs$tau
      #max_dim <- max( sapply(inputs$tau, function(tt) if(is.array(tt)) length(dim(tt)) else 1 ) )
      if(mixing) newDim <- c(inputs$nObs,
                             apollo_inputs$apollo_draws$interNDraws,
                             apollo_inputs$apollo_draws$intraNDraws)
      #if(max_dim==3 && intraNDraws==0) stop("Some thresholds for Ordered Logit for model component \"", inputs$componentName, "\" are given as arrays despite the model not using intra-individual mixing!")
      for(k in 1:length(inputs$tau)){ # expand if neccesary
        tt <- inputs$tau[[k]]
        if(is.function(tt)){ environment(tt) <- new.env(hash=TRUE, parent=parent.frame()); tt = tt() }
        if(is.vector(tt)){
          if(!(length(tt) %in% c(1, inputs$nObs))) stop("Threshold ", k, " for Ordered Logit for model component \"", inputs$componentName,"\" is given as a vector and should have one entry per observation in the data!")
          #if(length(tt)==1) inputs$tau[[k]] = rep(tt, inputs$nObs)  
          #if(max_dim   ==2) inputs$tau[[k]] = matrix(tt, nrow=newDim[1], ncol=newDim[2], byrow=FALSE)
          #if(max_dim   ==3) inputs$tau[[k]] = array(tt, dim=newDim)
        } else if(is.array(tt)){
          if(!mixing) stop("Model uses no mixing, so thresholds for Ordered Logit for model component \"",inputs$componentName,"\" should not be given as matrices or arrays!")
          isCube <- length(dim(tt))==3
          if(isCube && apollo_inputs$apollo_draws$intraNDraws==0) stop("Some thresholds for Ordered Logit for model component \"", inputs$componentName, "\" are given as arrays despite the model not using intra-individual mixing!")
          if(nrow(tt)!=inputs$nObs) stop("Threshold ",k," for Ordered Logit for model component \"",inputs$componentName,"\" is given as a matrix/array and should have one entry per observation in the data!")
          if(dim(tt)[2]!=newDim[2]) stop("Threshold ",k," for Ordered Logit for model component \"",inputs$componentName,"\" has a number of columns that is not equal to the number of inter-individual draws in the model!")
          #if(length(dim(tt))!=3) inputs$tau[[k]] = array(tt, dim=newDim)
          if(isCube && dim(tt)[3]!=newDim[3]) stop("Threshold ",k," for Ordered Logit for model component \"",inputs$componentName,"\" has a number of columns in the third dimension that is not equal to the number of intra-individual draws in the model!")
        } else stop("Some threshold has an unrecognised format for model component \"", inputs$componentName, "\".")
      }
    }
    
    if(is.null(inputs$coding)){
      inputs$coding <- 1:(length(inputs$tau)+1)  
      if(functionality=="validate"){
        txt <- paste0('No coding provided for Ordered Logit for model component "', inputs$componentName,
                      '", so assuming outcomeOrdered goes from 1 to ', max(inputs$coding))
        apollo_print(txt)
      }
    }
    
    ### Format checks
    # outcomeOrdered
    test <- is.vector(inputs$outcomeOrdered) && (length(inputs$outcomeOrdered)==inputs$nObs || length(inputs$outcomeOrdered)==1)
    if(!test) stop("The \"outcomeOrdered\" argument for model component \"",inputs$componentName,"\" needs to be a scalar or a vector with one entry per observation in the \"database\"")
    # V
    test <- is.numeric(inputs$V) || is.function(inputs$V)
    if(is.numeric(inputs$V)){
      test <- test && (is.vector(inputs$V) | is.array(inputs$V))
      lenV <- ifelse(is.array(inputs$V), dim(inputs$V)[1], length(inputs$V))
      test <- test && (lenV==inputs$nObs | lenV==1)
    } 
    if(!test) stop("\"V\" for model component \"",inputs$componentName,"\" needs to be a function, a scalar, or a vector/matrix/cube with one row per observation in the \"database\"")  
    # rows
    test <- is.vector(inputs$rows) && ( (is.logical(inputs$rows) && length(inputs$rows)==inputs$nObs) || (length(inputs$rows)==1 && inputs$rows=="all") )
    if(!test) stop("The \"rows\" argument for model component \"",inputs$componentName,"\" needs to be \"all\" or a vector of boolean statements with one entry per observation in the \"database\"")
    # coding
    test <- is.null(inputs$coding) || is.vector(inputs$coding)
    if(!test) stop("Argument 'coding', if provided for model component \"",inputs$componentName,"\", must be a vector.")
    
    # Expand rows if necessary, and update nObs
    if(length(inputs$rows)==1 && inputs$rows=="all") inputs$rows <- rep(TRUE, inputs$nObs)
    inputs$nObs <- sum(inputs$rows)
    # Filter rows, except for V
    if(any(!inputs$rows)){
      inputs$outcomeOrdered <- apollo_keepRows(inputs$outcomeOrdered, inputs$rows)
      inputs$tau <- lapply(inputs$tau, apollo_keepRows, inputs$rows)
    }
    
    # Apply coding
    map <- stats::setNames(1:length(inputs$coding), inputs$coding)
    inputs$outcomeOrdered <- map[as.character(inputs$outcomeOrdered)]
    
    # Add extreme thresholds NOT NECESSARY ANYMORE
    # inputs$tau <- c(-Inf, inputs$tau, Inf)
    
    # Create list of dummies for outcomeOrdered
    inputs$Y <- lapply(as.list(1:length(inputs$coding)), "==", inputs$outcomeOrdered)
    
  }
  
  #### MDCEV, MDCNEV ####
  if(modelType %in% c("mdcev", "mdcnev")){
    # Check for mandatory inputs
    mandatory <- c("alternatives", "continuousChoice", "V", "alpha", "gamma", "cost", "budget")
    if(modelType=="mdcev") mandatory <- c(mandatory, "sigma")
    if(modelType=="mdcnev") mandatory <- c(mandatory, "mdcnevNests", "mdcnevStructure")
    for(i in mandatory) if(!(i %in% names(inputs))) stop('The inputs list for model component "', inputs$componentName, 
                                                         '" needs to include an object called "', i,'"!')
    # Check for optional inputs
    if(is.null(inputs[["minConsumption"]])) inputs[["minConsumption"]] = NA
    if(is.null(inputs[["rows"]])          ) inputs[["rows"]]           = "all"
    if(is.null(inputs[["nRep"]])          ) inputs[["nRep"]]           = 100
    if(is.null(inputs[["outside"]])){
      if("outside" %in% inputs$alternatives) inputs$outside <- "outside" else inputs$outside <- NA
    } else if(!(inputs$outside %in% inputs$alternatives)) stop('Name provided for outside good for model component "',
                                                               inputs$componentName,'" does not correspond to an alternative!')
    if(modelType=="mdcnev"){
      test <- !is.null(inputs[["sigma"]]) && !(length(inputs$sigma)==1 && inputs$sigma==1) && !apollo_inputs$silent
      if(test) apollo_print(paste0('Setting "sigma" set to 1 in model component ', inputs$componentName,
                                   ' to ensure identifiability.'))
      inputs$sigma <- 1
    }
    if(is.null(inputs[['avail']])){
      inputs[['avail']]=1
      if(!apollo_inputs$silent && functionality=='validate') apollo_print('Setting "avail" is missing, so full availability is assumed.')
    }
    
    # continuousChoice
    test <- is.list(inputs$continuousChoice) && all(sapply(inputs$continuousChoice, is.numeric))
    if(!test ) stop('"continuousChoice" for model component "', inputs$componentName, 
                    '" should be a list of numerical vectors, each with as many rows as observations.')
    
    # V
    test <- is.list(inputs$V) && all(sapply(inputs$V, is.numeric))
    if(!test) stop('V for model component "', inputs$componentName, '" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.')
    
    # Useful variables
    inputs$nObs <- max(sapply(inputs$continuousChoice, length))
    inputs$nAlt <- length(inputs$V)
    
    # alternatives
    test <- is.vector(inputs$alternatives) && is.character(inputs$alternatives)
    if(!test) stop('"alternatives" for model component "', inputs$componentName, '" must be a character vector.')
    if(length(inputs$alternatives)!=inputs$nAlt) stop('The length of "alternatives" for model component "', 
                                                       inputs$componentName, '" does not match with the length of "V".')
    
    # avail
    inputs$avail_set <- FALSE
    if(length(inputs$avail)==1 && inputs$avail==1) { inputs$avail=rep(1, inputs$nObs); inputs$avail_set <- TRUE}
    test <- is.list(inputs$avail) && all(sapply(inputs$avail, function(x) is.logical(x) | is.numeric(x)))
    test <- test && all(sapply(inputs$avail, function(x) length(x) %in% c(1,inputs$nObs)))
    if(!test) stop('"avail" for model component "', inputs$componentName,
                   '" should be a list of scalars or vectors (each with as many elements as observations), with values 0 or 1.')
    lenA <- sapply(inputs$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
    test <- all(lenA==inputs$nObs | lenA==1)
    if(!test) stop('All entries in "avail" for model component "', inputs$componentName,
                   '" need to be a scalar or a vector with one entry per observation in the "database"')
    
    # alpha
    test <- is.list(inputs$alpha) && all(sapply(inputs$alpha, is.numeric)) 
    if(!test) stop('"alpha" for model component "', inputs$componentName,
                   '" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.')
    
    # gamma
    test <- is.list(inputs$gamma) && all(sapply(inputs$gamma, is.numeric))
    if(!test) stop('"gamma" for model component "', inputs$componentName, 
                   '" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.')
    
    # cost
    test <- is.list(inputs$cost) && all(sapply(inputs$cost, is.numeric))
    if(!test) stop('"cost" for model component "', inputs$componentName, 
                   '" should be a list of numerical scalar or vectors (each with as many rows as observations), with positive values.')
    
    # budget
    test <- is.numeric(inputs$budget) && ( length(inputs$budget) %in% c(1,inputs$nObs) )
    if(!test) stop('"budget" for model component "', inputs$componentName, 
                   '" should be a scalar or a vector with as many elements as observations.')
    
    # rows
    test <- is.vector(inputs$rows)
    test <- test && ( (is.logical(inputs$rows) && length(inputs$rows)==inputs$nObs) || (length(inputs$rows)==1 && inputs$rows=="all") )
    if(!test) stop('The "rows" argument for model component "', inputs$componentName,
                   '" needs to be "all" or a vector of logical statements with one entry per observation in the "database"')
    # minConsumption
    if(!anyNA(inputs$minConsumption)){
      test <- is.list(inputs$minConsumption) && length(inputs$minConsumption)!=inputs$nAlt
      test <- test && all(sapply(inputs$minConsumption, is.numeric)) && all(sapply(inputs$minConsumption, length) %in% c(1,inputs$nObs))
      if(!test) stop('"minConsumption" for model component "', inputs$componentName, 
                     '", if included, must be a list of numeric scalars or vectors (each with as many elements as observations).')
      # set minConsumption to NA if they are all zero
      if(Reduce('+',lapply(inputs$minConsumption, sum))==0) inputs$minConsumption=NA
    }
    if(anyNA(inputs$minConsumption)) inputs$minX <- FALSE else inputs$minX <- TRUE
    
    ## Add gamma outside, if missing
    test <- !anyNA(inputs$outside) && is.null(inputs$gamma[[inputs$outside]])
    if(test) inputs$gamma[[inputs$outside]] <- 1
    # Put the outsidegood first, if it exists
    if(!anyNA(inputs$outside)) inputs$alternatives=c(inputs$outside, 
                                                     inputs$alternatives[inputs$alternatives!=inputs$outside])
    # Sort lists based on the order of 'alternatives'
    #if(any(inputs$alternatives != names(inputs$V)               )) inputs$V                <- inputs$V[inputs$alternatives]
    if(any(inputs$alternatives != names(inputs$avail)           )) inputs$avail            <- inputs$avail[inputs$alternatives]
    #if(any(inputs$alternatives != names(inputs$alpha)           )) inputs$alpha            <- inputs$alpha[inputs$alternatives]
    #if(any(inputs$alternatives != names(inputs$gamma)           )) inputs$gamma            <- inputs$gamma[inputs$alternatives]
    if(any(inputs$alternatives != names(inputs$continuousChoice))) inputs$continuousChoice <- inputs$continuousChoice[inputs$alternatives]
    if(any(inputs$alternatives != names(inputs$cost)            )) inputs$cost             <- inputs$cost[inputs$alternatives]
    test <- !anyNA(inputs$minConsumption) && any(inputs$alternatives != names(inputs$minConsumption))
    if(test) inputs$minConsumption <- inputs$minConsumption[inputs$alternatives]
    # If there is an outside good, rename it to "outside"
    if(!anyNA(inputs$outside)){
      inputs$alternatives            = c("outside", inputs$alternatives[2:length(inputs$alternatives)])
      #names(inputs$V    )            = inputs$alternatives
      names(inputs$avail)            = inputs$alternatives
      #names(inputs$alpha)            = inputs$alternatives
      #names(inputs$gamma)            = inputs$alternatives
      names(inputs$cost )            = inputs$alternatives
      names(inputs$continuousChoice) = inputs$alternatives
      if(!anyNA(inputs$minConsumption)) names(inputs$minConsumption) = inputs$alternatives 
    }
    
    # Expand "rows"
    if(length(inputs$rows)==1 && inputs$rows=="all") inputs$rows <- rep(TRUE, inputs$nObs)
    # Remove excluded rows from cost, avail, continuousChoice, budget and minConsumption
    if(any(!inputs$rows)){
      inputs$cost             <- lapply(inputs$cost            , apollo_keepRows, r=inputs$rows)
      inputs$avail            <- lapply(inputs$avail           , apollo_keepRows, r=inputs$rows)
      inputs$continuousChoice <- lapply(inputs$continuousChoice, apollo_keepRows, r=inputs$rows)
      inputs$budget           <- apollo_keepRows(inputs$budget, r=inputs$rows)
      if(!anyNA(inputs$minConsumption)) inputs$minConsumption <- lapply(inputs$minConsumption, apollo_keepRows, r=inputs$rows)
    }
    
    # For unavailable alternatives: Set consumption to zero and cost to 1 (unless scalar).
    # No need no change cost if scalar, as an invalid number will fail for available rows too.
    for(j in 1:inputs$nAlt) if(any(!inputs$avail[[j]])){
      aSca <- length(inputs$avail[[j]])==1
      xSca <- length(inputs$continuousChoice[[j]])==1
      cSca <- length(inputs$cost[[j]])==1
      if(aSca){
        inputs$continuousChoice[[j]][!inputs$avail[[j]]] <- 0
        inputs$cost[[j]][!inputs$avail[[j]]] <- 1
      } else {
        if(xSca) inputs$continuousChoice[[j]] <- ifelse(inputs$avail[[j]], inputs$continuousChoice[[j]], 0) else inputs$continuousChoice[[j]][!inputs$avail[[j]]] <- 0
        if(!cSca) inputs$cost[[j]][!inputs$avail[[j]]] <- 1
      }
    }
    
    # Additional useful variables
    inputs$discrete_choice <- lapply(inputs$continuousChoice, ">", 0)
    inputs$totalChosen     <- Reduce("+",inputs$discrete_choice)
    inputs$chosenUnavail   <- mapply(function(m,a) m & !a, inputs$discrete_choice, inputs$avail, SIMPLIFY=FALSE)
    inputs$chosenUnavail   <- Reduce("|", inputs$choseUnavail)
    inputs$hasOutside      <- "outside" %in% inputs$alternatives
    inputs$altnames        <- inputs$alternatives
    if(modelType=="mdcnev"){
      inputs$nObs  <- nrow(apollo_inputs$database)
      if(any(!inputs$rows)) inputs$nObs <- sum(inputs$rows)
      inputs$nAlt <- length(inputs$V)
      inputs$nNests<- length(inputs$mdcnevNests)
      inputs$q = list() # stores how many different products where purchased in each nest
      for(s in 1:inputs$nNests){
        alts   <- which(as.vector(inputs$mdcnevStructure[s,])>0)
        inputs$q[[s]] <- Reduce("+", inputs$discrete_choice[alts])
      }; rm(alts)
      inputs$term4base=function(t, ars, qrs){
        sumx=1
        if(qrs==0|qrs==1) sumx = 1
        if(qrs >1 & ars==1) sumx = 1
        if (ars==2) sumx = ((1-t)/t)*qrs*(qrs-1)/2
        if (qrs==3 & ars==3){
          sumx = (2*(1-t)/t)+1
          sumx = sumx * ((1-t)/t)}
        if (qrs==4 & ars==3){
          sumx = ((3*(1-t)/t)+1)*((2*(1-t)/t)+(1*(1-t)/t))
          sumx = sumx + (((2*(1-t)/t)+1)*((1-t)/t))}
        if (qrs==4 & ars==4) sumx = (((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t))
        if (qrs==5 & ars==3){
          sumx = ((4*(1-t)/t)+1) * ((3*(1-t)/t)+(2*(1-t)/t)+(1*(1-t)/t))
          sumx = sumx + ((3*(1-t)/t)+1) * ((2*(1-t)/t)+(1*(1-t)/t))
          sumx = sumx + (((2*(1-t)/t)+1)*(1*(1-t)/t))}
        if (qrs==5 & ars==4){
          sumx = ((4*(1-t)/t)+2)*((3*(1-t)/t)+1) * ((2*(1-t)/t)+(1*(1-t)/t))
          sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)
          sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)}
        if (qrs==5 & ars==5) sumx = ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)
        if (qrs==6 & ars==3){
          sumx = ((5*(1-t)/t)+1)*(10*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==6 & ars==4){
          sumx = ((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==6 & ars==5){
          sumx = ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==6 & ars==6) sumx = ((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        if (qrs==7 & ars==3){
          sumx = ((6*(1-t)/t)+1)*(15*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+1)*(10*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==7 & ars==4){
          sumx = ((6*(1-t)/t)+2)*((5*(1-t)/t)+1)*(10*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(3*(1-t)/t)}
        if (qrs==7 & ars==5){
          sumx = ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==7 & ars==6){
          sumx = ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((6*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
          sumx = sumx + ((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
        if (qrs==7 & ars==7) sumx = ((6*(1-t)/t)+5)*((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        return(sumx)
      }
      inputs$term4final=function(mdcnevNests, r_current_combo, q_person){
        out <- 1
        for(s in 1:length(r_current_combo)) if(q_person[s]>0) out <- out * inputs$term4base(mdcnevNests[[s]], r_current_combo[s], q_person[s])
        return(out)
      }
      inputs$term5=function(r_current_combo, q_person){
        out <- 0
        for(s in 1:length(r_current_combo)) if(q_person[s]>0) out <- out + q_person[s]- r_current_combo[s] + 1
        out <- factorial(out - 1)
        return(out)
      }
      
    }
    
    
  }
  
  return(inputs)
  
  
}
