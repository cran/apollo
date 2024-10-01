#' Pre-process input for common models
#' return
#' @param inputs List of settings
#' @param modelType Character. Type of model, e.g. "mnl", "nl", "cnl", etc.
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
#' @param apollo_inputs List of main inputs to the model estimation process. See \link{apollo_validateInputs}.
#' @return The returned object depends on the value of argument operation
#' @export
apollo_validate <- function(inputs, modelType, functionality, apollo_inputs){
  modelType <- tolower(modelType)
  
  #### MNL, NL, CNL, DFT, EL ####
  if(modelType %in% c("mnl","fmnl","nl","cnl","dft", "el")){
    
    # Check there are no repeated alternatives names
    if(length(unique(inputs$altnames))!=length(inputs$altnames)) stop('SPECIFICATION ISSUE - Names of alternatives must be unique. Check definition of "alternatives".')
    
    # Check that no alternative is called "chosen" or "choice"
    test <- any(tolower(inputs$altnames) %in% c("chosen", "choice"))
    if(test) stop("SYNTAX ISSUE - Alternatives cannot be called 'chosen' or 'choice'!")
    
    # Checks specific to CNL
    if(modelType=="cnl"){
      if("root" %in% names(inputs$cnlNests)) stop('SPECIFICATION ISSUE - The root should not be included in argument cnlNests for model component "', inputs$componentName,'".')
      test <- all(inputs$cnlNests!=0)
      if(!test) stop("SPECIFICATION ISSUE - Structural parameters (lambda) cannot be zero for model component \"",inputs$componentName,"\"!")      
      test <- is.matrix(inputs$cnlStructure) && nrow(inputs$cnlStructure)==length(inputs$nestnames) && ncol(inputs$cnlStructure)==inputs$nAlt
      if(!test) stop('SPECIFICATION ISSUE - Argument "cnlStructure" for model component "', inputs$componentName, ' must be a matrix with one row per nest and one column per alternative.')
      #test <- 0.999<colSums(inputs$cnlStructure) && colSums(inputs$cnlStructure)<1.001
      test <- all((0.999<colSums(inputs$cnlStructure)) & (colSums(inputs$cnlStructure)<1.001))
      if(!test) stop("SPECIFICATION ISSUE - Allocation parameters (alpha) for some alternatives sum to values different than 1 for model component \"",inputs$componentName,"\"!")
      test <- all((0<=inputs$cnlStructure) & (inputs$cnlStructure<=1))
      if(!test) stop("SPECIFICATION ISSUE - Some allocation parameters (alpha) are outside the 0-1 range for model component \"",inputs$componentName,"\"!")
    }
    
    # Check that there are at least two or three alternatives
    minAlts <- 2; if(modelType%in%c("el","nl","cnl")) minAlts <- 3
    if(inputs$nAlt<minAlts) stop("SPECIFICATION ISSUE - Model component \"",inputs$componentName,"\"  requires at least ", minAlts, " alternatives")
    
    # Check that choice vector is not empty
    if(modelType!="fmnl"){
      if(modelType!="el") if(length(inputs$choiceVar)==0) stop("SPECIFICATION ISSUE - Choice vector is empty for model component \"",inputs$componentName,"\"")
      if(modelType=="el") for(i in 1:length(inputs$choiceVars)){
        if(length(inputs$choiceVars[[i]])==0) stop('SPECIFICATION ISSUE - Choice vector is empty for stage ',i,' in model component "',inputs$componentName,'"')
      } 
    } else {
      for(i in 1:length(inputs$choiceShares)){
        if(length(inputs$choiceShares[[i]])==0) stop('SPECIFICATION ISSUE - Choice shares vector is empty for alternative ',inputs$altnames[i],' in model component "',inputs$componentName,'"')}
    }
    if(inputs$nObs==0) stop("INPUT ISSUE - No data for model component \"",inputs$componentName,"\"")
    
    # Check that shares are all between 0 and 1, and sum to 1 for FMNL
    if(modelType=="fmnl"){
      test <- any(sapply(inputs$choiceShares,anyNA))
      if(any(test)) stop('INPUT ISSUE - Choice shares for some observations for model component "', inputs$componentName,
                         '" are NA!')
      test <- any(sapply(inputs$choiceShares,">",1))
      if(any(test)) stop('INPUT ISSUE - Choice shares for some observations for model component "', inputs$componentName,
                         '" are greater than 1!')
      test <- any(sapply(inputs$choiceShares,"<",0))
      if(any(test)) stop('INPUT ISSUE - Choice shares for some observations for model component "', inputs$componentName,
                         '" are less than 0!')
      totalChoice <- Reduce("+",inputs$choiceShares)
      test <- (abs(totalChoice - 1) > 0.001)
      if(any(test)) stop('INPUT ISSUE - Choice shares for some observations for model component "', inputs$componentName,
             '" do not sum to 1!')
    }

    # Check V and avail elements are named correctly
    if(modelType!="dft" && !all(inputs$altnames %in% names(inputs$V))) stop("SPECIFICATION ISSUE - The names of the alternatives for model component \"",inputs$componentName,"\" do not match those in \"utilities\".")
    if(modelType!="el") if(!all(inputs$altnames %in% names(inputs$avail))) stop("SPECIFICATION ISSUE - The names of the alternatives for model component \"",inputs$componentName,"\" do not match those in \"avail\".")
    if(modelType=="el") for(s in 1:inputs$stages) if(!all(inputs$altnames %in% names(inputs$avail[[s]]))) stop('SPECIFICATION ISSUE - The names of the alternatives for model component "',inputs$componentName,'" do not match those in "avail" (in stage ',s,').')
    
    # Check that there are no values in the choice column for undefined alternatives
    if(modelType!="fmnl"){
      if(modelType!="el"){
        inputs$choiceLabs <- unique(inputs$choiceVar)
        if(!all(inputs$choiceLabs %in% inputs$altcodes)) stop("SPECIFICATION ISSUE - The data contains values for \"choiceVar\" for model component \"",inputs$componentName,"\" that are not included in \"alternatives\".")
      } else {
        choiceLabs <- unique(unlist(inputs$choiceVars))
        if(!all(choiceLabs %in% c(inputs$altcodes,-1))) stop("SPECIFICATION ISSUE - The data contains values for \"choiceVar\" for model component \"",inputs$componentName,"\" that are not included in \"alternatives\".")
      }
    }
    
    # Checks specific for Exploded Logit (EL)
    if(modelType=="el"){
      # check that all availabilities are either 0 or 1
      for(i in 1:length(inputs$avail)) if( !all(unlist(inputs$avail[[i]]) %in% 0:1) ) stop("INPUT ISSUE - Some availability values for model component \"",inputs$componentName,"\" are not 0 or 1.")
      # check that at least 2 alternatives are available in at least one observation
      for(i in 1:length(inputs$avail)) if(max(Reduce('+',inputs$avail[[i]]))==1) stop("INPUT ISSUE - Only one alternative is available for each observation for model component \"",inputs$componentName,"!")
      # check that nothing unavailable is chosen
      for(s in 1:inputs$stages) for(j in 1:inputs$nAlt){
        tmp <- !inputs$avail[[s]][[j]] & inputs$Y[[s]][[j]]
        if(any(tmp) && any(!inputs$rows)) tmp <- apollo_insertRows(tmp, inputs$rows, FALSE)
        tmp <- paste0(which(tmp), collapse=',')
        if(nchar(tmp)>0 && !apollo_inputs$silent) stop(paste0('INPUT ISSUE - Alternative "', inputs$altnames[j], 
                                                              '" is chosen in row(s) ', tmp, ' in stage ', s, ' of model ',
                                                              'component "', inputs$componentName, 
                                                              '", despite not being available.'))
      }
    } else {
      # check that all availabilities are either 0 or 1
      for(i in 1:length(inputs$avail)) if( !all(unique(inputs$avail[[i]]) %in% 0:1) ) stop("SPECIFICATION ISSUE - Some availability values for model component \"",inputs$componentName,"\" are not 0 or 1.")
      # check that at least 2 alternatives are available in at least one observation
      if(max(Reduce('+',inputs$avail))==1) stop("INPUT ISSUE - Only one alternative is available for each observation for model component \"",inputs$componentName,"!")
      # check that nothing unavailable is chosen
      if(modelType!="fmnl"){
        for(j in 1:inputs$nAlt) if(any(inputs$choiceVar==inputs$altcodes[j] & inputs$avail[[j]]==0)){
          #stop("INPUT ISSUE - The data contains cases where alternative ",
          #     inputs$altnames[j]," is chosen for model component \"",
          #     inputs$componentName,"\" despite being listed as unavailable\n")
          txt <-  paste0('The data contains cases where alternative ', inputs$altnames[j], 
                         ' is chosen for model component "',inputs$componentName, '" despite being', 
                         ' listed as unavailable. This will cause the chosen probability to be', 
                         ' zero, and potentially lead to an invalid LL.')
          apollo_print(txt, type="w")
        } 
      } else {
        for(j in 1:inputs$nAlt) if(any(inputs$choiceShares[[j]]>0 & inputs$avail[[j]]==0)){
          #stop("The data contains cases where alternative ", inputs$altnames[j], 
          #     " is chosen for model component \"", inputs$componentName, 
          #     "\" despite being listed as unavailable\n")
          txt <-  paste0('The data contains cases where alternative ', inputs$altnames[j], 
                         ' is chosen for model component "',inputs$componentName, '" despite being', 
                         ' listed as unavailable. This will cause the chosen probability to be', 
                         ' zero, and potentially lead to an invalid LL.')
          apollo_print(txt, type="w")
        }
      }
    }
    
    # Check that no available alternative has utility = NA
    # Requires setting non available alternatives utility to 0 first
    if(modelType=="el") inputs$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V, inputs$avail[[1]], SIMPLIFY=FALSE)
    if(!(modelType %in% c("el", "dft"))) inputs$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V, inputs$avail, SIMPLIFY=FALSE)
    if(modelType!="dft" && !all(sapply(inputs$V, function(v) all(is.finite(v))))) stop('CALCULATION ISSUE - Some utilities for model component "',
                                                                                       inputs$componentName, 
                                                                                       '" contain values that are not finite numbers!')
    
    if(modelType=='nl'){
      allElements <- c("root", unlist(inputs$nlStructure))
      if(is.null(inputs$nlStructure[["root"]])) stop("SPECIFICATION ISSUE - Tree structure for model component \"",inputs$componentName,"\" is missing an element called root!")
      #test <- all(inputs$nlNests!=0)
      test <- all(sapply(inputs$nlNests, function(x) is.numeric(x) && all(x!=0)))
      if(!test) stop("SPECIFICATION ISSUE - Structural parameters (lambda) cannot be zero for model component \"",inputs$componentName,"\"!")      
      #if(inputs$nlNests[["root"]]!=1) stop("SPECIFICATION ISSUE - The root lambda parameter for model component \"",inputs$componentName,"\" should be equal to 1.")
      if( !all(inputs$altnames %in% allElements) ) stop("SPECIFICATION ISSUE - All alternatives must be included in the tree structure for model component \"",inputs$componentName,"\".")
      if( !all(inputs$nestnames %in% allElements) ) stop("SPECIFICATION ISSUE - All nests must be included in the tree structure for model component \"",inputs$componentName,"\".")
      if( (length(inputs$nestnames)+length(inputs$altnames))!=length(allElements) ) stop("SPECIFICATION ISSUE - Tree structure for model component \"",inputs$componentName,"\" is inconsistent. Each element must appear only once.")
      if( !all(names(inputs$nlNests) %in% names(inputs$nlStructure)) | !all(names(inputs$nlStructure) %in% names(inputs$nlNests)) ) stop("SPECIFICATION ISSUE - All nests in argument 'nlNests' for model component \"",inputs$componentName,"\" should be in 'nlStructure', and vice versa (including 'root').")
      combined_elements="root"
      for(j in 1:length(inputs$nlStructure)) combined_elements=c(combined_elements, inputs$nlStructure[[j]])
      for(j in 1:length(inputs$altnames)){
        if(sum(inputs$nestnames==inputs$altnames[j])) stop("SPECIFICATION ISSUE - A nest for model component \"",inputs$componentName,"\" cannot have the same name as an alternative!")
        if(sum(combined_elements==inputs$altnames[j])!=1) stop("SPECIFICATION ISSUE - An alternative for model component \"",inputs$componentName,"\" needs to appear exactly once in a tree!")
      }
      for(j in 1:length(inputs$nlStructure)) if(sum(inputs$nestnames==names(inputs$nlStructure)[j])!=1) stop("SPECIFICATION ISSUE - A defined nest for model component \"",
                                                                                                             inputs$componentName,"\" needs to appear exactly once in a tree!")
      for(j in 1:length(inputs$nestnames)){
        if(sum(inputs$altnames==inputs$nestnames[j])) stop("SPECIFICATION ISSUE - A nest for model component \"",inputs$componentName,"\" cannot have the same name as an alternative!")
        if(sum(combined_elements==inputs$nestnames[j])!=1) stop("SPECIFICATION ISSUE - A defined nest for model component \"",inputs$componentName,"\" needs to appear exactly once in a tree!")
      }
    } # end of NL specific validation
    
  }
  
  #### classAlloc #### 
  if(modelType=='classalloc'){
    # Check there are at least 2 classes
    if(inputs$nAlt<2) stop("SPECIFICATION ISSUE - Model component \"",inputs$componentName,"\"  requires at least 2 classes")
    # Check availabilities are only 0 or 1
    if(!all(sapply(inputs$avail, function(a) all(unique(a) %in% 0:1)))) stop("SPECIFICATION ISSUE - Some availability values for model component \"",inputs$componentName,"\" are not 0 or 1.")
    # Check all available alternatives have finite V
    inputs$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V, inputs$avail, SIMPLIFY=FALSE)
    test     <- all(sapply(inputs$V, function(v) all(is.finite(v))))
    if(!test) stop('SPECIFICATION ISSUE - Some utilities for model component "',inputs$componentName, '" contain values that are not finite numbers!')
    # Check that the name in classes, V and avail match
    ### retrieve class names (as the vector classes could be character or a named vector)
    if(is.character(inputs$classes)){
      inputs$classNames=inputs$classes
    }else{
      inputs$classNames=names(inputs$classes)
    }
    test <- length(inputs$V)==length(inputs$classes) && length(inputs$V)==length(inputs$avail)
    test <- test && all(names(inputs$V) %in% inputs$classNames) && all(names(inputs$V) %in% names(inputs$avail))
    if(!test) stop("SPECIFICATION ISSUE - Some of the names of the classes for model component '", inputs$componentName, "' do not match!")
  }
  
  #### NormD ####
  if(modelType=="normd"){
    if(is.vector(inputs$xNormal)) xlength=length(inputs$xNormal)
    if(is.array(inputs$xNormal)) xlength=dim(inputs$xNormal)[1]
    if(!is.vector(inputs$outcomeNormal)) stop("INPUT ISSUE - Dependent variable for model component \"",inputs$componentName,"\" needs to be one-dimensional!")
    if(xlength!=1 && xlength!=length(inputs$outcomeNormal)) stop("INPUT ISSUE - Incompatible dimensions for dependent and explanatory variables for model component \"",inputs$componentName,"\"!")
    if(!all(is.finite(inputs$xNormal))) stop('INPUT ISSUE - Some values inside xNormal are not finite for model component "', inputs$componentName, '"')
  }
  
  #### NormD ####
  if(modelType=="tobit"){
    if(is.vector(inputs$xTobit)) xlength=length(inputs$xTobit)
    if(is.array(inputs$xTobit)) xlength=dim(inputs$xTobit)[1]
    if(!is.vector(inputs$outcomeTobit)) stop("INPUT ISSUE - Dependent variable for model component \"",inputs$componentName,"\" needs to be one-dimensional!")
    if(xlength!=1 && xlength!=length(inputs$outcomeTobit)) stop("INPUT ISSUE - Incompatible dimensions for dependent and explanatory variables for model component \"",inputs$componentName,"\"!")
    if(!all(is.finite(inputs$xTobit))) stop('INPUT ISSUE - Some values inside xTobit are not finite for model component "', inputs$componentName, '"')
    if(length(inputs$lowerLimit)!=1) stop("INPUT ISSUE - Lower limit for model component \"",inputs$componentName,"\" needs to be a scalar!")
    if(length(inputs$upperLimit)!=1) stop("INPUT ISSUE - Upper limit for model component \"",inputs$componentName,"\" needs to be a scalar!") 
    if(any(inputs$outcomeTobit<inputs$lowerLimit)) stop("INPUT ISSUE - Dependent variable for model component \"",inputs$componentName,"\" below lower limit for some observations!")
    if(any(inputs$outcomeTobit>inputs$upperLimit)) stop("INPUT ISSUE - Dependent variable for model component \"",inputs$componentName,"\" above upper limit for some observations!")
  }

  #### OL, OP ####
  if(modelType %in% c("ol", "op")){
    # validation of coding is done in apollo_preprocess, as coding is applied in there
    if( (length(inputs$tau)+1) != length(inputs$coding) ) stop("SPECIFICATION ISSUE - Threshold vector length +1 does not match number of elements in argument 'coding' for model component \"",inputs$componentName,"\".")
    if(length(inputs$tau)>1) for(s in 1:(length(inputs$tau)-1)) if(any(inputs$tau[[s+1]]<=inputs$tau[[s]])) stop('SPECIFICATION ISSUE - Tresholds for model component "', inputs$componentName, '" are not strictly increasing at starting values!')
    if(!all(is.finite(inputs$V))) stop('CALCULATION ISSUE - Some values inside V are not finite for model component "', inputs$componentName, '"')
  }
  
  #### MDCEV, MDCNEV ####
  if(modelType %in% c("mdcev", "mdcnev")){
    # Check names of alternatives are unique
    if(length(inputs$alternatives)!=length(unique(inputs$alternatives))) stop('SPECIFICATION ISSUE - Alternatives names must be unique. Check definition of "alternatives".')
    # Check that sigma is not random (actually, it could be, but it leads to weird results)
    if(!is.vector(inputs$sigma)) stop("SPECIFICATION ISSUE - Sigma for model component \"", inputs$componentName,"\" should not be random")
    if(!(length(inputs$sigma) %in% c(1,inputs$nObs))) stop("SPECIFICATION ISSUE - Sigma for model component \"",inputs$componentName,"\" should be either a scalar or a vector with as many elements as observations")
    # Check there are at least two alternatives
    if(inputs$nAlt<2) stop("SPECIFICATION ISSUE - Model component \"",inputs$componentName,"\" requires at least two products")
    # Check that choice vector is not empty
    if(inputs$nObs==0) stop("INPUT ISSUE - No data for model component \"",inputs$componentName,"\"")
    # Check that MDCNEV has outside good
    if(!inputs$hasOutside && modelType=="mdcnev") stop("SPECIFICATION ISSUE - The MDCNEV structured used for model component \"",inputs$componentName,"\" requires an \"outside\" good, i.e. an alternative for which consumption is always non-zero!")
    # Check that first product is outside good (if one is defined)
    if(inputs$hasOutside && inputs$alternatives[1]!="outside") stop("SPECIFICATION ISSUE - First product for model component \"",inputs$componentName,"\" must be called \"outside\"!")
    # Check labels
    namesinside=names(inputs$V)[names(inputs$V)!="outside"]
    if(!all(inputs$alternatives %in% names(inputs$V))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"utilities\"!")
    if(!all(inputs$alternatives %in% names(inputs$alpha))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"alpha\"!")
    if(!all(namesinside %in% names(inputs$gamma))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"gamma\"!")
    if(!all(inputs$alternatives %in% names(inputs$continuousChoice))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"continuousChoice\"!")
    if(!all(inputs$alternatives %in% names(inputs$cost))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"cost\"!")
    if(!all(inputs$alternatives %in% names(inputs$avail))) stop("SPECIFICATION ISSUE - Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"avail\"!")
    # check that nothing unavailable is chosen
    for(j in 1:inputs$nAlt) if( any( inputs$discrete_choice[[j]] & inputs$avail[[j]]==0 ) ) stop("INPUT ISSUE - Product", inputs$alternatives[j], "chosen despite being listed as unavailable for model component \"", inputs$componentName,"\"!")
    # check that outside good is always chosen
    if(inputs$hasOutside){
      txt <- paste0("INPUT ISSUE - Outside good ", ifelse(inputs$outside!="outside", paste0("(",inputs$outside,") "), ""),
                    "for model component \"", inputs$componentName,"\" should always be chosen!")
      if(any(inputs$continuousChoice[["outside"]]<=0)) stop(txt)
    }
    # check that all costs are positive
    if( sum(sapply(inputs$cost, function(x) sum(x<=0))) > 0 ) stop("INPUT ISSUE - Costs for model component \"",inputs$componentName,
                                                                   "\" must be strictly positive for all products!")
    # check consumption is non negative for all products
    for(i in 1:length(inputs$continuousChoice)) if(any(inputs$continuousChoice[[i]]<0)){
      stop( paste0("INPUT ISSUE - Consumption values of alternative", inputs$alternatives[i], 
                   "for model component \"",inputs$componentName,"\" must be non-negative!") )
    }
    # Check budget>0
    if(any(inputs$budget<=0)) stop("INPUT ISSUE - Budget for model component \"",inputs$componentName,
                                   "\" for some rows in data is less than or equal to zero!")
    # check that full budget is consumed in each row, nothing more, nothing less
    expenditure <- Reduce("+", mapply("*", inputs$continuousChoice, inputs$cost, SIMPLIFY=FALSE))
    test <- which(abs(expenditure/inputs$budget - 1) > 0.001)
    if(length(test)>0){
      df <- data.frame(ID     = apollo_inputs$database[test, apollo_inputs$apollo_control$indivID],
                       budget = if(length(inputs$budget)==1) inputs$budget else inputs$budget[test], 
                       expend = expenditure[test],
                       `%diff`= round((expenditure[test]/inputs$budget[test] - 1)*100, 2),
                       check.names = FALSE)
      rownames(df) <- test
      df <- df[order(abs(df[,4]), decreasing=TRUE),]
      print(df)
      stop('INPUT ISSUE - Expenditure for some observations for model component "', inputs$componentName,
           '" is either less or more than budget!')
    }
    # turn scalar availabilities into vectors
    #for(i in 1:length(inputs$avail)) if(length(inputs$avail[[i]])==1) inputs$avail[[i]] <- rep(inputs$avail[[i]], inputs$nObs)
    # check that all availabilities are either 0 or 1
    for(i in 1:inputs$nAlt) if( !all(unique(inputs$avail[[i]]) %in% 0:1) ) stop("INPUT ISSUE - Some availability values are not 0 or 1 for model component \"",
                                                                                inputs$componentName,"\".")
    # check that availability of outside is always 1
    if(inputs$hasOutside && any(!inputs$avail[[1]])) stop('INPUT ISSUE - Outside good is not available for some observations for model component "',
                                                          inputs$componentName,'". It should always be available.')
    # check that if minimum consumption exists, it has the same names as alternatives, and that no consumptions are less than minConsumption if alternative is available
    if(inputs$minX){
      if(!all(inputs$alternatives %in% names(inputs$minConsumption))) stop("SYNTAX ISSUE - Labels in \"alternatives\" for model component \"",
                                                                           inputs$componentName,"\" do not match those in \"minConsumption\"!")
      for(i in 1:inputs$nAlt){
        test <- any(inputs$continuousChoice[[i]][inputs$avail[[i]]] < inputs$minConsumption[[i]][inputs$avail[[i]]])
        if(test) stop( paste0("INPUT ISSUE - Consumption of alternative ", inputs$alternatives[i], " for model component \"",
                              inputs$componentName,"\" is smaller than its listed minConsumption") )
      }
    }
    
    # Checks specific to MDCNEV
    if(modelType=="mdcnev"){
      if(nrow(inputs$mdcnevStructure)!=length(inputs$mdcnevNests)) stop("SPECIFICATION ISSUE - Tree structure needs one row per nest!")
      if(ncol(inputs$mdcnevStructure)!=inputs$nAlt) stop("SPECIFICATION ISSUE - Tree structure needs one column per alternative!")
      if(any(colSums(inputs$mdcnevStructure)!=1)) stop("ESPECIFICATION ISSUE - ach alternative must be allocated to one nest only!")
      if(!all(inputs$mdcnevStructure %in% 0:1)) stop("SPECIFICATION ISSUE - All values in mdcnevStructure must be either 0 or 1!")
    }
    
  }
  
  return(TRUE)
}
