#' Pre-process input for multiple models
#' return
#' @param inputs List of settings
#' @param modelType Character. Type of model, e.g. "mnl", "nl", "cnl", etc.
#' @param functionality Character. Either "estimate","prediction","validate","zero_LL","conditionals","output","raw", or "preprocess". Only used for validation, it does not influence the return values.
#' @param apollo_inputs List of main inputs to the model estimation process. See \link{apollo_validateInputs}.
#' @return The returned object depends on the value of argument operation
#' @export
apollo_validate <- function(inputs, modelType, functionality, apollo_inputs){
  modelType <- tolower(modelType)
  
  #### MNL, NL, CNL, DFT, EL ####
  if(modelType %in% c("mnl","nl","cnl","dft", "el")){
    
    # Checks specific to CNL
    if(modelType=="cnl"){
      if("root" %in% names(inputs$cnlNests)) stop('The root should not be included in argument cnlNests for model component "', inputs$componentName,'".')
      test <- is.matrix(inputs$cnlStructure) && nrow(inputs$cnlStructure)==length(inputs$nestnames) && ncol(inputs$cnlStructure)==inputs$nAlt
      if(!test) stop('Argument "cnlStructure" for model component "', inputs$componentName, ' must be a matrix with one row per nest and one column per alternative.')
      test <- 0.999<colSums(inputs$cnlStructure) && colSums(inputs$cnlStructure)<1.001
      if(!test) stop("Allocation parameters (alpha) for some alternatives sum to values different than 1 for model component \"",inputs$componentName,"\"!")
    }
    
    # Check that there are at least two or three alternatives
    minAlts <- 2; if(modelType%in%c("el","nl","cnl")) minAlts <- 3
    if(inputs$nAlt<minAlts) stop("Model component \"",inputs$componentName,"\"  requires at least ", minAlts, " alternatives")
    
    # Check that choice vector is not empty
    if(modelType!="el") if(length(inputs$choiceVar)==0) stop("Choice vector is empty for model component \"",inputs$componentName,"\"")
    if(modelType=="el") for(i in 1:length(inputs$choiceVars)){
      if(length(inputs$choiceVars[[i]])==0) stop('Choice vector is empty for stage ',i,' in model component "',inputs$componentName,'"')
    } 
    if(inputs$nObs==0) stop("No data for model component \"",inputs$componentName,"\"")
    
    # Check V and avail elements are named correctly
    if(modelType!="dft" && !all(inputs$altnames %in% names(inputs$V))) stop("The names of the alternatives for model component \"",inputs$componentName,"\" do not match those in \"V\".")
    if(modelType!="el") if(!all(inputs$altnames %in% names(inputs$avail))) stop("The names of the alternatives for model component \"",inputs$componentName,"\" do not match those in \"avail\".")
    if(modelType=="el") for(s in 1:inputs$stages) if(!all(inputs$altnames %in% names(inputs$avail[[s]]))) stop('The names of the alternatives for model component "',inputs$componentName,'" do not match those in "avail" (in stage ',s,').')
    
    # Check that there are no values in the choice column for undefined alternatives
    if(modelType!="el"){
      inputs$choiceLabs <- unique(inputs$choiceVar)
      if(!all(inputs$choiceLabs %in% inputs$altcodes)) stop("The data contains values for \"choiceVar\" for model component \"",inputs$componentName,"\" that are not included in \"alternatives\".")
    } else {
      choiceLabs <- unique(unlist(inputs$choiceVars))
      if(!all(choiceLabs %in% c(inputs$altcodes,-1))) stop("The data contains values for \"choiceVar\" for model component \"",inputs$componentName,"\" that are not included in \"alternatives\".")
    }
    
    # Checks specific for Exploded Logit (EL)
    if(modelType=="el"){
      # check that all availabilities are either 0 or 1
      for(i in 1:length(inputs$avail)) if( !all(unlist(inputs$avail[[i]]) %in% 0:1) ) stop("Some availability values for model component \"",inputs$componentName,"\" are not 0 or 1.")
      # check that nothing unavailable is chosen
      for(s in 1:inputs$stages) for(j in 1:inputs$nAlt){
        tmp <- !inputs$avail[[s]][[j]] & inputs$Y[[s]][[j]]
        if(any(tmp) && any(!inputs$rows)) tmp <- apollo_insertRows(tmp, inputs$rows, FALSE)
        tmp <- paste0(which(tmp), collapse=',')
        if(nchar(tmp)>0 && !apollo_inputs$silent) stop(paste0('Alternative "', inputs$altnames[j], 
                                                              '" is chosen in row(s) ', tmp, ' in stage ', s, ' of model ',
                                                              'component "', inputs$componentName, 
                                                              '", despite not being available.'))
      }
    } else {
      # check that all availabilities are either 0 or 1
      for(i in 1:length(inputs$avail)) if( !all(unique(inputs$avail[[i]]) %in% 0:1) ) stop("Some availability values for model component \"",inputs$componentName,"\" are not 0 or 1.")
      # check that nothing unavailable is chosen
      for(j in 1:inputs$nAlt) if(any(inputs$choiceVar==inputs$altcodes[j] & inputs$avail[[j]]==0)) stop("The data contains cases where alternative ",
                                                                                                        inputs$altnames[j]," is chosen for model component \"",
                                                                                                        inputs$componentName,"\" despite being listed as unavailable\n")
    }
    
    # Check that no available alternative has utility = NA
    # Requires setting non available alternatives utility to 0 first
    if(modelType=="el") inputs$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V, inputs$avail[[1]], SIMPLIFY=FALSE)
    if(!(modelType %in% c("el", "dft"))) inputs$V <- mapply(function(v,a) apollo_setRows(v, !a, 0), inputs$V, inputs$avail, SIMPLIFY=FALSE)
    if(modelType!="dft" && !all(sapply(inputs$V, function(v) all(is.finite(v))))) stop('Some utilities for model component "',
                                                                                       inputs$componentName, 
                                                                                       '" contain values that are not finite numbers!')
    
    if(modelType=='nl'){
      allElements <- c("root", unlist(inputs$nlStructure))
      if(is.null(inputs$nlStructure[["root"]])) stop("Tree structure for model component \"",inputs$componentName,"\" is missing an element called root!")
      if(inputs$nlNests["root"]!=1) stop("The root lambda parameter for model component \"",inputs$componentName,"\" should be equal to 1.")
      if( !all(inputs$altnames %in% allElements) ) stop("All alternatives must be included in the tree structure for model component \"",inputs$componentName,"\".")
      if( !all(inputs$nestnames %in% allElements) ) stop("All nests must be included in the tree structure for model component \"",inputs$componentName,"\".")
      if( (length(inputs$nestnames)+length(inputs$altnames))!=length(allElements) ) stop("Tree structure for model component \"",inputs$componentName,"\" is inconsistent. Each element must appear only once.")
      if( !all(names(inputs$nlNests) %in% names(inputs$nlStructure)) | !all(names(inputs$nlStructure) %in% names(inputs$nlNests)) ) stop("All nests in argument 'nlNests' for model component \"",inputs$componentName,"\" should be in 'nlStructure', and vice versa (including 'root').")
      combined_elements="root"
      for(j in 1:length(inputs$nlStructure)) combined_elements=c(combined_elements, inputs$nlStructure[[j]])
      for(j in 1:length(inputs$altnames)){
        if(sum(inputs$nestnames==inputs$altnames[j])) stop("A nest for model component \"",inputs$componentName,"\" cannot have the same name as an alternative!")
        if(sum(combined_elements==inputs$altnames[j])!=1) stop("An alternative for model component \"",inputs$componentName,"\" needs to appear exactly once in a tree!")
      }
      for(j in 1:length(inputs$nlStructure)) if(sum(inputs$nestnames==names(inputs$nlStructure)[j])!=1) stop("A defined nest for model component \"",
                                                                                                             inputs$componentName,"\" needs to appear exactly once in a tree!")
      for(j in 1:length(inputs$nestnames)){
        if(sum(inputs$altnames==inputs$nestnames[j])) stop("A nest for model component \"",inputs$componentName,"\" cannot have the same name as an alternative!")
        if(sum(combined_elements==inputs$nestnames[j])!=1) stop("A defined nest for model component \"",inputs$componentName,"\" needs to appear exactly once in a tree!")
      }
    } # end of NL specific validation
    
  }
  
  #### NormD ####
  if(modelType=="normd"){
    if(is.vector(inputs$xNormal)) xlength=length(inputs$xNormal)
    if(is.array(inputs$xNormal)) xlength=dim(inputs$xNormal)[1]
    if(!is.vector(inputs$outcomeNormal)) stop("Dependent variable for model component \"",inputs$componentName,"\" needs to be one-dimensional!")
    if(xlength!=1 && xlength!=length(inputs$outcomeNormal)) stop("Incompatible dimensions for dependent and explanatory variables for model component \"",inputs$componentName,"\"!")
    if(!all(is.finite(inputs$xNormal))) stop('Some values inside xNormal are not finite for model component "', inputs$componentName, '"')
  }
  
  #### OL, OP ####
  if(modelType %in% c("ol", "op")){
    values_present = unique(inputs$outcomeOrdered)
    if(!(all(values_present %in% inputs$coding ))) stop("Some levels in 'outcomeOrdered' do not exist in 'coding' for model component \"",inputs$componentName,"\" !")
    if(!(all(inputs$coding %in% values_present ))) stop("Some levels in 'coding' do not exist in 'outcomeOrdered' for model component \"",inputs$componentName,"\"!")
    if( (length(inputs$tau)+1) != length(inputs$coding) ) stop("Threshold vector length +1 does not match number of elements in argument 'coding' for model component \"",inputs$componentName,"\".")
    if(!all(is.finite(inputs$V))) stop('Some values inside V are not finite for model component "', inputs$componentName, '"')
  }
  
  #### MDCEV ####
  if(modelType %in% c("mdcev", "mdcnev")){ 
    # Check that sigma is not random (actually, it could be, but it leads to weird results)
    if(!is.vector(inputs$sigma)) stop("Sigma for model component \"", inputs$componentName,"\" should not be random")
    if(!(length(inputs$sigma) %in% c(1,inputs$nObs))) stop("Sigma for model component \"",inputs$componentName,"\" should be either a scalar or a vector with as many elements as observations")
    # Check there are at least two alternatives
    if(inputs$nAlt<2) stop("Model component \"",inputs$componentName,"\" requires at least two products")
    # Check that choice vector is not empty
    if(inputs$nObs==0) stop("No data for model component \"",inputs$componentName,"\"")
    # Check that first product is outside good
    if(inputs$hasOutside && inputs$alternatives[1]!="outside") stop("First product for model component \"",inputs$componentName,"\" must be called \"outside\"!")
    # Check labels
    namesinside=names(inputs$V)[names(inputs$V)!="outside"]
    if(!all(inputs$alternatives %in% names(inputs$V))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"V\"!")
    if(!all(inputs$alternatives %in% names(inputs$alpha))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"alpha\"!")
    if(!all(namesinside %in% names(inputs$gamma))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"gamma\"!")
    if(!all(inputs$alternatives %in% names(inputs$continuousChoice))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"continuousChoice\"!")
    if(!all(inputs$alternatives %in% names(inputs$cost))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"cost\"!")
    if(!all(inputs$alternatives %in% names(inputs$avail))) stop("Labels in \"alternatives\" for model component \"",inputs$componentName,"\" do not match those in \"avail\"!")
    # check that nothing unavailable is chosen
    for(j in 1:inputs$nAlt) if( any( inputs$discrete_choice[[j]] & inputs$avail[[j]]==0 ) ) stop("Product", inputs$alternatives[j], "chosen despite being listed as unavailable for model component \"", inputs$componentName,"\"!")
    # check that outside good is always chosen
    if(inputs$hasOutside){
      txt <- paste0("Outside good ", ifelse(inputs$outside!="outside", paste0("(",inputs$outside,") "), ""),
                    "for model component \"", inputs$componentName,"\" should always be chosen!")
      if(any(inputs$continuousChoice[["outside"]]<=0)) stop(txt)
    }
    # check that all costs are positive
    if( sum(sapply(inputs$cost, function(x) sum(x<=0))) > 0 ) stop("Costs for model component \"",inputs$componentName,
                                                                   "\" must be strictly positive for all products!")
    # check consumption is non negative for all products
    for(i in 1:length(inputs$continuousChoice)) if(any(inputs$continuousChoice[[i]]<0)){
      stop( paste0("Consumption values of alternative", inputs$alternatives[i], 
                   "for model component \"",inputs$componentName,"\" must be non-negative!") )
    }
    # Check budget>0
    if(any(inputs$budget<=0)) stop("Budget for model component \"",inputs$componentName,
                                   "\" for some rows in data is less than or equal to zero!")
    # check that full budget is consumed in each row, nothing more, nothing less
    expenditure <- Reduce("+", mapply("*", inputs$continuousChoice, inputs$cost, SIMPLIFY=FALSE))
    if(any(abs(expenditure-inputs$budget)>10^-10)) stop("Expenditure for some observations for model component \"",
                                                        inputs$componentName,"\" is either less or more than budget!")
    # turn scalar availabilities into vectors
    #for(i in 1:length(inputs$avail)) if(length(inputs$avail[[i]])==1) inputs$avail[[i]] <- rep(inputs$avail[[i]], inputs$nObs)
    # check that all availabilities are either 0 or 1
    for(i in 1:inputs$nAlt) if( !all(unique(inputs$avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1 for model component \"",
                                                                         inputs$componentName,"\".")
    # check that availability of outside is always 1
    if(inputs$hasOutside && any(!inputs$avail[[1]])) stop('Outside good is not available for some observations for model component "',
                                                   inputs$componentName,'". It should always be available.')
    # check that if minimum consumption exists, it has the same names as alternatives, and that no consumptions are less than minConsumption if alternative is available
    if(inputs$minX){
      if(!all(inputs$alternatives %in% names(inputs$minConsumption))) stop("Labels in \"alternatives\" for model component \"",
                                                                           inputs$componentName,"\" do not match those in \"minConsumption\"!")
      for(i in 1:inputs$nAlt){
        test <- any(inputs$continuousChoice[[i]][inputs$avail[[i]]] < inputs$minConsumption[[i]][inputs$avail[[i]]])
        if(test) stop( paste0("Consumption of alternative ", inputs$alternatives[i], " for model component \"",
                              inputs$componentName,"\" is smaller than its listed minConsumption") )
      }
    }
    
    # Checks specific to MDCNEV
    if(modelType=="mdcnev"){
      if(nrow(inputs$mdcnevStructure)!=length(inputs$mdcnevNests)) stop("Tree structure needs one row per nest!")
      if(ncol(inputs$mdcnevStructure)!=inputs$nAlt) stop("Tree structure needs one column per alternative!")
      if(any(colSums(inputs$mdcnevStructure)!=1)) stop("Each alternative must be allocated to one nest only!")
    }
    
  }
  
  return(TRUE)
}
