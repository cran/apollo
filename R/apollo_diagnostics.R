#' Pre-process input for common models
#' return
#' @param inputs List of settings
#' @param modelType Character. Type of model, e.g. "mnl", "nl", "cnl", etc.
#' @param apollo_inputs List of main inputs to the model estimation process. See \link{apollo_validateInputs}.
#' @param data Boolean. TRUE for printing report related to dependent and independent variables. FALSE for not printing it. Default is TRUE.
#' @param param Boolean. TRUE for printing report related to estimated parameters (e.g. model structure). FALSE for not printing it. Default is TRUE.
#' @return (invisibly) TRUE if no error happend during execution.
#' @export
apollo_diagnostics <- function(inputs, modelType, apollo_inputs, data=TRUE, param=TRUE){
  modelType <- tolower(modelType)
  
  #### MNL, NL, CNL, DFT ####
  if(modelType %in% c("mnl","nl","cnl","dft")){
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
      if(any(choicematrix[4,]==1)) apollo_print("Some alternatives are always chosen when available!", type="w")
      #if(inputs$avail_set) apollo_print("Availability not provided (or some elements are NA). Full availability assumed.", type="w")
      apollo_print("\n")
      apollo_print(paste0('Overview of choices for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      print(round(choicematrix,2))
    }
    
    if(modelType=="cnl" & param){
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
        apollo_print(paste0('Structure for ', toupper(modelType), ' model component ', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
        apollo_print(out_tree)
        for(a in 1:(ncol(out_tree)-1)){
          if(sum(out_tree[,a])!=1) apollo_print(paste0("Allocation parameters for alternative \'",inputs$altnames[a],"\' do not sum to 1. You may want to impose a constraint using a logistic transform."), type="i")
          if(any(out_tree[,a]<0)) apollo_print(paste0("Some allocation parameters for alternative \'",inputs$altnames[a],"\' are negative. You may want to impose a constraint using a logistic transform."), type="i")
          if(any(out_tree[,a]>1)) apollo_print(paste0("Some allocation parameters for alternative \'",inputs$altnames[a],"\' are larger than 1. You may want to impose a constraint using a logistic transform."), type="i")
        }
        if(any(out_tree[,ncol(out_tree)]<0)) apollo_print("Some lambda parameters are negative. You may want to impose a constraint or reconsider your model structure.", type="i")
        if(any(out_tree[,ncol(out_tree)]>1)) apollo_print("Some lambda parameters are larger than 1. You may want to impose a constraint or reconsider your model structure.", type="i")
      }
    }
    
    if(modelType=='nl' & param){
      if(!apollo_inputs$silent & data) apollo_print('\n') # 
      if(!apollo_inputs$silent){
        # WARNING for automatic setting of root nesting parameter
        if(inputs$root_set) apollo_print("Root lambda parameter set to 1.", type="i")
        # Identifying nest's parents
        nestAbove <- unique(lapply(inputs$ancestors, '[', -1))
        nestAbove <- setNames(sapply(nestAbove, function(x) if(length(x)==1) return('Inf') else x[2]) ,
                              sapply(nestAbove, '[', 1))
        # Printing graphical representation of the tree, using recursive function
        apollo_print(paste0('Nesting structure for ', toupper(modelType), ' model component ', 
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
              cat("\n",space,"-Nest: ", nlStructure[[component]][j], " (",round(l,4), ")", sep="")
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
    } # end of NL special checks
    
    if(modelType=="dft" & !apollo_inputs$silent & data){
      txt <- 'Not all of the attributes given in "attrValues" are named in "attrScalings" or "attrWeights". These will consequently be ignored.'
      if(inputs$warn1) apollo_print(txt, type="i")
      txt <- 'Not all of the alternatives given in "altStart" are named in "alternatives". These will consequently be ignored.'
      if(inputs$warn2) apollo_print(txt, type="i")
      txt <- 'A list was not supplied for "altStart". Starting values for all alternatives will be set to zero.'
      if(inputs$warn3) apollo_print(txt, type="i")
    }
    
  } 
  
  #### EL ####
  if(modelType=="el"){
    ### changes 28 July
    # Initialise summary table of availabilities and market share
    #choicematrix <- array(0, dim=c(4, inputs$nAlt+1, inputs$stages), 
    #                      dimnames=list(c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available"),
    #                                    c(inputs$altnames, "No choice"), paste("stage", 1:inputs$stages)) )
    choicematrix <- array(0, dim=c(4, inputs$nAlt, inputs$stages), 
                          dimnames=list(c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available"),
                                        c(inputs$altnames), paste("stage", 1:inputs$stages)) )
    # Calculate summary table for each stage and print it
    for(s in 1:inputs$stages){
      for(a in 1:inputs$nAlt){
        choicematrix[1,a,s] <- ifelse(length(inputs$avail[[s]][[a]])==1 && inputs$avail[[s]][[a]]==1, 
                                      inputs$nObs, sum(inputs$avail[[s]][[a]]) )
        choicematrix[2,a,s] <- sum(inputs$Y[[s]][[a]])
        choicematrix[3,a,s] <- choicematrix[2,a,s]/inputs$nObs*100
        choicematrix[4,a,s] <- choicematrix[2,a,s]/choicematrix[1,a,s]*100
        if(!is.finite(choicematrix[4,a,s])) choicematrix[4,a,s] <- 0
      }
      if(!apollo_inputs$silent & data){
        apollo_print("\n")
        apollo_print(paste0('Overview of choices for ', toupper(modelType), ' model component', 
                            ifelse(inputs$componentName=='model', '', inputs$componentName), ', stage ', s,':'))
        print(round(choicematrix[,,s],2))
      }
    }
    
    if(!apollo_inputs$silent & data) for(a in 1:inputs$nAlt){
      if(sum(choicematrix[4,a,])==0) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'), type="w")
      if(choicematrix[4,a,1]==1) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always chosen when available in model component "', inputs$componentName, '".'), type="w")
    }
    #if(inputs$avail_set==TRUE & !apollo_inputs$silent & data) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'), type="i")
  }
  
  #### FMNL ####
  if(modelType=="fmnl"){
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
      apollo_print(paste0('Overview of choices for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      print(round(choicematrix,2))
      
      # Print warnings
      for(a in 1:inputs$nAlt){
        if(choicematrix[2,a]==0) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'), type="w")
        if(choicematrix[4,a]==1) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always given the full choice share when available in model component "', inputs$componentName, '".'), type="w")
      }
      #if(inputs$avail_set==TRUE & !apollo_inputs$silent) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'), type="i")
    }
  }
  
  #### NormD ####
  if(modelType=="normd"){
    if(!apollo_inputs$silent & data){
      apollo_print('\n')
      apollo_print(paste0('Summary statistics for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      tmp <- t(as.matrix(summary(inputs$outcomeNormal))); rownames(tmp) <- ""
      print(round(tmp,4))
      rm(tmp)
    }
  }
  
  #### OL, OP ####
  if(modelType %in% c("ol", "op")){
    choicematrix <- t(as.matrix(table(inputs$outcomeOrdered)))
    choicematrix <- rbind(choicematrix, choicematrix[1,]/inputs$nObs*100)
    rownames(choicematrix) <- c("Times chosen", "Percentage chosen overall")
    if(!apollo_inputs$silent & data){
      apollo_print("\n")
      apollo_print(paste0('Overview of choices for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      print(round(choicematrix,2))
    }
  }
  
  #### MDCEV, MDCNEV ####
  if(modelType %in% c("mdcev", "mdcnev")){
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
      apollo_print(paste0('Overview of choices for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      print(round(choicematrix,2))
      
      # Print warnings
      for(a in 1:inputs$nAlt){
        if(choicematrix[2,a]==0) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is never chosen in model component "', inputs$componentName, '".'), type="w")
        if(choicematrix[2,a]==choicematrix[1,a] && inputs$altnames[a]!=inputs$outside) apollo_print(paste0('Alternative "', inputs$altnames[a], '" is always chosen when available in model component "', inputs$componentName, '".'), type="w")
      }
      #if(inputs$avail_set==TRUE & !apollo_inputs$silent) apollo_print(paste0('Availability not provided (or some elements are NA) for model component ', inputs$componentName,'. Full availability assumed.'))
    }
    
    if(modelType=="mdcnev" && !apollo_inputs$silent & param){
      if(data) apollo_print('\n')
      apollo_print(paste0('Nest structure for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      colnames(inputs$mdcnevStructure) <- inputs$alternatives
      rownames(inputs$mdcnevStructure) <- names(inputs$mdcnevNests)
      maxL <- max(nchar(names(inputs$mdcnevNests)))
      for(n in rownames(inputs$mdcnevStructure)) apollo_print(paste0(
        n, paste0(rep('.', maxL - nchar(n)), collapse=''),' (', round(inputs$mdcnevNests[[n]], 2), '): ', 
        paste0(inputs$alternatives[inputs$mdcnevStructure[n,]>0], collapse=', ')
      ))
      #apollo_print(inputs$mdcnevStructure)
    }
  }
  
  #### LC ####
  if(modelType=='lc'){
    class_summary = matrix(0, nrow=length(inputs$classProb), ncol=1, dimnames=list(names(inputs$inClassProb), "Mean prob."))
    if(is.null(rownames(class_summary))) rownames(class_summary) <- paste0("Class_", 1:length(inputs$classProb))
    for(cc in 1:length(inputs$classProb)) class_summary[cc,1] <- mean(inputs$classProb[[cc]])
    if(!apollo_inputs$silent & param){
      apollo_print('\n')
      apollo_print(paste0('Summary of class allocation for ', toupper(modelType), ' model component ', 
                          ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
      apollo_print(class_summary)
    }
  }
  
  return(invisible(TRUE))
}
