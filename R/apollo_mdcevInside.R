#' Calculates MDCEV likelihoods without an outside good.
#' 
#' Calculates the likelihood of a Multiple Discrete Continuous Extreme Value (MDCEV) model without an outside good.
#' 
#' @param V Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#' @param alternatives Character vector. Names of alternatives, elements must match the names in list 'V'.
#' @param alpha Named list. Alpha parameters for each alternative. As many elements as alternatives.
#' @param gamma Named list. Gamma parameters for each alternative. As many elements as alternatives.
#' @param sigma Numeric scalar. Scale parameter of the model extreme value type I error.
#' @param cost Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#' @param avail Named list. Availabilities of alternatives, one element per alternative. Names of elements must match those in argument 'alternatives'. Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation. If all alternatives are always available, then user can just omit this argument.
#' @param continuousChoice Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#' @param budget Numeric vector. Budget for each observation.
#' @param minConsumption Named list of scalars or numeric vectors. Minimum consumption of the alternatives, if consumed. As many elements as alternatives. Names must match those in \code{alternatives}.
#' @param rows Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate" Used for model estimation.
#'                        \item "prediction" Used for model predictions.
#'                        \item "validate" Used for validating input.
#'                        \item "zero_LL" Used for calculating null likelihood.
#'                        \item "conditionals" Used for calculating conditionals.
#'                        \item "output" Used for preparing output after model estimation.
#'                        \item "raw" Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": A matrix with one row per observation, and means and s.d. of predicted consumptions.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": Not applicable.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @importFrom stats setNames
apollo_mdcevInside <- function(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice,budget,functionality,minConsumption=NA, rows="all"){
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Store useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    avail_set <- FALSE
    
    # check rows statement
    if(length(rows)!=nObs & !(length(rows)==1 && rows=="all")) stop("The argument \"rows\" needs to either be \"all\" or a boolean vector of length equal to the number of the rows in the data!")
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[1]][!rows] = budget[!rows] else {
        continuousChoice[[1]][!rows] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }
    
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), alternatives)
    }
    
    discrete_choice=list()
    k=1
    while(k<= nAlts){
      discrete_choice[[alternatives[k]]]=(continuousChoice[[k]]>0)
      k=k+1
    }
    
    # set minConsumption to NA if they are all zero
    if(!anyNA(minConsumption)) if(Reduce('+',lapply(minConsumption, sum))==0) minConsumption=NA
    
    apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
                                error=function(e) return(list(noValidation=FALSE, noDiagnostics=FALSE)) )
    
    if(apollo_control$noValidation==FALSE){
      # Check that sigma is not random (actually, it could be, but it leads to weird results)
      if(!is.vector(sigma)) stop("Sigma should not be random")
      if(length(sigma)!=1 && length(sigma)!=nObs) stop("Sigma should be either a scalar or a vector with as many elements as observations")
      
      # Check there are at least two alternatives
      if(nAlts<2) stop("MDCEV requires at least two products")
      
      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")  
      
      
      # Check labels
      if(!all(alternatives %in% names(V))) stop("Labels in \"alternatives\" do not match those in \"V\"!")
      if(!all(alternatives %in% names(alpha))) stop("Labels in \"alternatives\" do not match those in \"alpha\"!")
      if(!all(alternatives %in% names(gamma))) stop("Labels in \"alternatives\" do not match those in \"gamma\"!")
      if(!all(alternatives %in% names(continuousChoice))) stop("Labels in \"alternatives\" do not match those in \"continuousChoice\"!")
      if(!all(alternatives %in% names(cost))) stop("Labels in \"alternatives\" do not match those in \"cost\"!")
      if(!all(alternatives %in% names(avail))) stop("Labels in \"alternatives\" do not match those in \"avail\"!")
      
      # check that nothing unavailable is chosen
      chosenunavail=0
      j=1
      while(j<= length(V)){
        if(sum(discrete_choice[[j]]*(avail[[j]]==0))) chosenunavail=1
        j=j+1}
      if(chosenunavail==1) stop("Some product(s) chosen despite being listed as unavailable!")
      
      # checks for individual products
      budget_check=0*budget
      
      k=1
      while(k<= length(V)){
        # check that all costs are positive
        if(sum(cost[[k]]==0)>0) stop("Need strictly positive costs for all products!")
        
        # check that no negative consumptions for any products
        if(sum(continuousChoice[[k]]<0)>0) stop("Negative consumption for some products for some observations!")
        budget_check=budget_check+continuousChoice[[k]]*cost[[k]]
        k=k+1
      }
      
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")
      
      # check that if minimum consumption exists, it has the same names as alternatives, and that no consumptions are less than minConsumption if alternative is available
      if(!anyNA(minConsumption)){
        if(!all(alternatives %in% names(minConsumption))) stop("Labels in \"alternatives\" do not match those in \"minConsumption\"!")
        # first ensure order is correct
        if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
        if(any(alternatives != names(minConsumption))) minConsumption <- minConsumption[alternatives]
        consumption_below_min_flag=0
        j=1
        while(j<=length(alternatives)){
          consumption_below_min_flag=consumption_below_min_flag+sum((continuousChoice[[j]]>0)*(continuousChoice[[j]]<minConsumption[[j]]))
          j=j+1
        }
        if(consumption_below_min_flag) stop("\nSome consumptions are below the lower limits listed in \"minConsumption\"!")
      }
      
      # confirm all checks are passed, print output and return TRUE
      #cat("\nAll checks passed for MDCEV model component\n")
    }
    
    if(apollo_control$noDiagnostics==FALSE){
      choicematrix = matrix(0,nrow=4,ncol=length(V))
      choicematrix[1,] = unlist(lapply(avail, function(x) sum(x[rows])))
      #choicematrix[1,] = colSums(matrix(unlist(avail), ncol = length(avail))[rows,])
      j=1
      while(j<= length(V)){
        choicematrix[2,j]=sum(discrete_choice[[j]][rows]==1)
        choicematrix[3,j]=sum(continuousChoice[[j]][rows])/choicematrix[1,j]
        choicematrix[4,j]=sum(continuousChoice[[j]][rows])/choicematrix[2,j]
        j=j+1
      }
      choicematrix <- apply(choicematrix, MARGIN=2, function(x) {x[!is.finite(x)] <- 0; return(x)})
      rownames(choicematrix) = c("Times available","Observations in which chosen","Average consumption when available","Average consumption when chosen")
      colnames(choicematrix) = names(V)
      
      #cat('Overview of choices for MDEV model component:\n')
      #print(round(choicematrix,2))
      #cat("\n")
      #if(any(choicematrix[2,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      #if(any(choicematrix[2,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                           "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for MDCEV model component:", content, apolloLog)
    }
    
    return(TRUE)
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    nObs  <- length(continuousChoice[[1]])
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[1]][!rows] = budget[!rows] else {
        continuousChoice[[1]][!rows] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }
    P <- rep(NA,nObs)
    P[!rows] <- 1
    return(P)
  }
  
  
  
  
  # ############################################### #
  #### functionality="estimate/conditionals/raw" ####
  # ################################################# #
  
  if(functionality%in%c("estimate","conditionals","raw")){
    # Store useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    # check rows statement
    if(length(rows)!=nObs & !(length(rows)==1 && rows=="all")) stop("The argument \"rows\" needs to either be \"all\" or a boolean vector of length equal to the number of the rows in the data!")
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[j]][1] = budget[!rows] else {
        continuousChoice[[j]][1] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), alternatives)
    }
    # Get discrete_choice
    discrete_choice=list()
    k=1
    while(k<= nAlts){
      discrete_choice[[alternatives[k]]]=(continuousChoice[[k]]>0)
      k=k+1
    }
    # set minConsumption to NA if they are all zero
    if(!anyNA(minConsumption)) if(Reduce('+',lapply(minConsumption, sum))==0) minConsumption=NA
    
    # Reorder V and avail to match alternatives order, if necessary
    if(any(alternatives != names(V))) V <- V[alternatives]
    if(any(alternatives != names(avail))) avail <- avail[alternatives]
    if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
    if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
    if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
    if(any(alternatives != names(discrete_choice))) discrete_choice <- discrete_choice[alternatives]
    if(any(alternatives != names(cost))) cost <- cost[alternatives]
    if(!anyNA(minConsumption)) if(any(alternatives != names(minConsumption))) minConsumption <- minConsumption[alternatives]
    
    # Compute V
    j=1      
    while(j<=length(V)){
      if(!anyNA(minConsumption)){
        tmp <- continuousChoice[[j]]-(continuousChoice[[j]]>=minConsumption[[j]])*minConsumption[[j]]
        V[[j]] = V[[j]] + avail[[j]]*((alpha[[j]]-1)*log((tmp/gamma[[j]]) + 1) - log(cost[[j]]))
      } else {
        V[[j]] = V[[j]] + avail[[j]]*((alpha[[j]]-1)*log((continuousChoice[[j]]/gamma[[j]]) + 1) - log(cost[[j]]))
      }
      j=j+1
    }
    
    # Get M (number of different chosen products)
    totalChosen=Reduce("+",discrete_choice)
    
    # first term
    term1=(1-totalChosen)*log(sigma)
    
    # compute log of fi for inside goods for use in second term, this has one column per inside good
    logfi=list()
    j=1   
    while(j<=length(V)){
      if(!anyNA(minConsumption)){
        tmp <- continuousChoice[[j]]-(continuousChoice[[j]]>=minConsumption[[j]])*minConsumption[[j]]
        logfi[[j]]=avail[[j]]*( log(1-alpha[[j]]) - log(tmp + gamma[[j]]) )     
      } else {
        logfi[[j]]=avail[[j]]*(log(1-alpha[[j]])-log(continuousChoice[[j]]+gamma[[j]]))  
      }
      j=j+1
    }
    
    # second term
    term2=0  
    j=1  
    while(j<=length(V)){
      term2=term2+avail[[j]]*(logfi[[j]]*discrete_choice[[j]])   
      j=j+1
    }
    
    # Third term
    term3 = 0  
    j=1       
    while(j<=length(V)){
      term3=term3+avail[[j]]*(cost[[j]]/exp(logfi[[j]]) * discrete_choice[[j]])   
      j=j+1}
    term3 = log(term3)
    
    # Fourth term
    term4_1 = 0
    term4_2 = 0
    j=1        
    while(j<=length(V)){
      term4_1 = term4_1+avail[[j]]*(V[[j]]/sigma * discrete_choice[[j]])
      term4_2 = term4_2+avail[[j]]*exp(V[[j]]/sigma)
      j=j+1
    }
    term4_2 = totalChosen * log(term4_2) # log of the above to the power of M (totalChosen)
    term4 =  term4_1 - term4_2 # log of fourth term is now the difference of the two parts
    
    # fifth term: log of factorial
    term5 = lfactorial(totalChosen-1)
    
    # probability is simply the exp of the sum of the logs of the individual terms
    P = exp(term1 + term2 + term3 + term4 + term5)
    
    if(is.vector(P)) P[!rows]   <- 1
    if(is.matrix(P)) P[!rows,]  <- 1
    if(is.array(P) && length(dim(P))==3) P[!rows,,] <- 1
    
    return(P)
  }
  
  
  # ############################## #
  #### functionality="output" ####
  # ############################## #
  
  if(functionality=="output"){
    # Calculate likelihood
    P <- apollo_mdcevInside(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice, budget, functionality="estimate", minConsumption, rows)
    
    # Useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    avail_set <- FALSE
    # Check 'rows' command
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[1]][!rows] = budget[!rows] else {
        continuousChoice[[1]][!rows] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), alternatives)
    }
    
    discrete_choice=list()
    for(k in 1:nAlts) discrete_choice[[alternatives[k]]] = (continuousChoice[[k]]>0)*1
    if(!anyNA(minConsumption)) if(Reduce('+',lapply(minConsumption, sum))==0) minConsumption=NA # set minConsumption to NA if they are all zero
    if(any(alternatives != names(V))) V <- V[alternatives] # Reorder V and avail to match alternatives order, if necessary
    if(any(alternatives != names(avail))) avail <- avail[alternatives]
    if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
    if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
    if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
    if(any(alternatives != names(cost))) cost <- cost[alternatives]
    
    # turn scalar availabilities into vectors
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
    
    availprint = colSums(matrix(unlist(avail), ncol = length(avail))[rows,])
    
    choicematrix = matrix(0,nrow=4,ncol=length(V))
    
    choicematrix[1,] = availprint
    j=1
    while(j<= length(V)){
      choicematrix[2,j]=sum(discrete_choice[[j]][rows]==1)
      choicematrix[3,j]=sum(continuousChoice[[j]][rows])/choicematrix[1,j]
      choicematrix[4,j]=sum(continuousChoice[[j]][rows])/choicematrix[2,j]
      j=j+1
    }
    choicematrix <- apply(choicematrix, MARGIN=2, function(x) {x[!is.finite(x)] <- 0; return(x)})
    rownames(choicematrix) = c("Times available","Observations in which chosen","Average consumption when available","Average consumption when chosen")
    colnames(choicematrix) = names(V)
    
    ## write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
    #apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE )$apollo_control,
    #                            error=function(e){
    #                              cat("apollo_mdcev could not retrieve apollo_control. No diagnostics in output.\n")
    #                              return(NA)
    #                            } )
    #if(!(length(apollo_control)==1 && is.na(apollo_control))){
    #  fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    #  fileName <- file.path(tempdir(),fileName)
    #  fileConn <- tryCatch( file(fileName, open="at"),
    #                        error=function(e){
    #                          cat('apollo_mdcev could not write diagnostics to temporary file. No diagnostics in output.\n')
    #                          return(NA)
    #                        })
    #  if(!anyNA(fileConn)){
    #    sink(fileConn)
    #    on.exit({if(sink.number()>0) sink(); close(fileConn)})
    #    if(apollo_control$noDiagnostics==FALSE){
    #      cat('\nOverview of choices for MDEV model component:\n')
    #      print(round(choicematrix,2))
    #      cat('\n')}
    #    if(sum(choicematrix[2,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
    #  }
    #}
    
    content <- list(round(choicematrix,2))
    if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for MDCEV model component:", content, apolloLog)
    
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #
  
  if(functionality=="prediction"){
    
    # Useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[1]][!rows] = budget[!rows] else {
        continuousChoice[[1]][!rows] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }
    discrete_choice=list()
    for(k in 1:nAlts) discrete_choice[[alternatives[k]]] = (continuousChoice[[k]]>0)*1
    
    # set minConsumption to NA if they are all zero
    if(!anyNA(minConsumption)) if(Reduce('+',lapply(minConsumption, sum))==0) minConsumption=NA
    
    # Reorder V and avail to match alternatives order, if necessary
    if(any(alternatives != names(V))) V <- V[alternatives]
    if(any(alternatives != names(avail))) avail <- avail[alternatives]
    if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
    if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
    if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
    if(any(alternatives != names(cost))) cost <- cost[alternatives]
    
    # first include a check to make sure that all alphas are equal across alts
    equality_check=1
    j=2
    while(j<= length(alpha)){
      if(alpha[[j]]!=alpha[[1]]) equality_check=0
      j=j+1
    }
    if(equality_check!=1) stop("\nMDCEV prediction only implemented for profile where alpha is generic across all products!")
    
    # include a check to make sure there are no minimum consumption limits
    if(!anyNA(minConsumption)) stop("\nMDCEV prediction only implemented for case without non-zero minimum consumptions!")
    
    forecasting=function(V,alternatives,alpha,gamma,sigma,cost,avail,budget,continuousChoice)
    {
      set.seed(99)
      
      ND=250
      nObs=length(continuousChoice[[1]])[1]
      consumption_total=array(0,dim=c(nObs,length(alternatives),ND))
      
      # create uniform draws for value draws
      # this needs to change still if we also have random alpha, gamma or V
      
      Ndraws=apollo_mlhs(ND*nObs,length(alternatives),1)
      draws_base=cbind(x=rep(seq(1,ND),nObs),Ndraws)
          
      k=1
      cat("\n0%")
      while(k<(ND+1))
      {
        draws=subset(draws_base,draws_base[,1]==k)
        draws=draws[,2:ncol(draws)]
        
        draws=sigma*(-log(-log(draws)))      
        
        V_trans=V
        
        
        j=1                                           
        while(j<=length(V)){
          V_trans[[j]]=exp(V_trans[[j]]+draws[,j])
          V_trans[[j]]=V_trans[[j]]/(cost[[j]])
          V_trans[[j]]=V_trans[[j]]*avail[[j]]
          j=j+1
        }
        
        for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
        for(i in 1:length(alpha)) if(length(alpha[[i]])==1) alpha[[i]] <- rep(alpha[[i]], nObs)
        for(i in 1:length(gamma)) if(length(gamma[[i]])==1) gamma[[i]] <- rep(gamma[[i]], nObs)
        for(i in 1:length(cost)) if(length(cost[[i]])==1) cost[[i]] <- rep(cost[[i]], nObs)
               
        # turn into matrices
        V_trans_use=Reduce("cbind",V_trans)
        alpha_use=Reduce("cbind",alpha)
        gamma_use=Reduce("cbind",gamma)
        avail_use=Reduce("cbind",avail)
        cost_use=Reduce("cbind",cost)
        
        consumption_total[,,k]=forecasting_subroutine(V_trans_use,alpha_use,gamma_use,avail_use,cost_use,budget,continuousChoice)
        
        # increment draws  
        if(k%%round(ND/10,0)==0) if(k==round(ND/2,0)) cat("50%") else if(k==ND) cat("100%") else cat(".")
        k=k+1
      }
      cat("\n")
      return(consumption_total)
    }
    
    forecasting_subroutine=function(V_trans,alpha,gamma,avail,cost,budget,continuousChoice)
    {
      nObs=length(continuousChoice[[1]])[1]
      N=length(V_trans[,1])
      NALT=ncol(V_trans)
      consumption=array(0,dim=c(N,NALT))
      i=1
      
      while(i<(N+1)){
        orderofV=rank(-V_trans[i,1:NALT])     
        M=1
        stopping=0
        while(stopping<1)
        {
          use=(orderofV<=M)   
          #step2
          lambda_1=((budget[i]+sum(cost[i,(1:ncol(cost))]*gamma[i,(1:ncol(gamma))]*use)))   
          lambda_21 = 0                                   
          lambda_22=sum(cost[i,(1:ncol(cost))]*gamma[i,(1:ncol(gamma))]*use*(V_trans[i,1:NALT])^(1/(1-alpha[i,1])))  
          lambda=(lambda_1/(lambda_21+lambda_22))^(alpha[i,1]-1)
          if((sum(lambda<V_trans[i,1:NALT]))<=(M))   
          {              
            #step3
            consumption_inside_1=(V_trans[i,1:NALT])^(1/(1-alpha[i,1:NALT]))*(budget[i]+sum(use*cost[i,(1:ncol(cost))]*gamma[i,(1:ncol(gamma))]))
            consumption_inside_2=lambda_21+lambda_22
            consumption[i,1:NALT]=use*((consumption_inside_1/consumption_inside_2)-1)*gamma[i,(1:ncol(gamma))]  
            stopping=1
          }
          else
          {
            #step4
            M=M+1  
            use=(orderofV<=M)   
            if(M==(sum(avail[i,1:NALT])))  
            {
              lambda_22=sum(cost[i,(1:ncol(cost))]*gamma[i,(1:ncol(gamma))]*use*(V_trans[i,1:NALT])^(1/(1-alpha[i,1])))  
              consumption_inside_1=(V_trans[i,1:NALT])^(1/(1-alpha[i,1:NALT]))*(budget[i]+sum(use*cost[i,(1:ncol(cost))]*gamma[i,(1:ncol(gamma))]))
              consumption_inside_2=lambda_21+lambda_22
              consumption[i,1:NALT]=use*((consumption_inside_1/consumption_inside_2)-1)*gamma[i,(1:ncol(gamma))]  
              stopping=1  
            }
          }
        }
        i=i+1
      } 
      return(consumption)
    }
    
    apollo_inputs <- tryCatch( get("apollo_inputs", parent.frame(), inherits=TRUE ),
                               error=function(e){
                                 cat("apollo_mdcev could not retrieve apollo_inputs.\n")
                                 cat(" Assuming no mixing.\n")
                                 return( list(apollo_control=list(mixing=FALSE)) )
                               } )
    apollo_control <- apollo_inputs$apollo_control
    if(apollo_control$mixing==FALSE){
      cat("\nNow producing forecasts from MDCEV model component.")
      cat("\nA matrix with one row per observation and the following columns will be returned:")
      cat("\n1. Means of predicted continuous consumptions (one column per product)")
      cat("\n2. Std err of predicted continuous consumptions (one column per product)")
      cat("\n3. Means of predicted discrete consumptions (one column per product)")
      cat("\n4. Std err of predicted discrete consumptions (one column per product)")
      cat("\n\nThis may take a while!")
      
      consumption=forecasting(V,alternatives,alpha,gamma,sigma,cost,avail,budget,continuousChoice)
      
      dimnames(consumption)[[2]]=alternatives
      
      predictions_continuous_mean=apply(consumption,c(1,2),mean)
      predictions_continuous_sd=apply(consumption,c(1,2),stats::sd)
      predictions_discrete=(consumption>0)
      predictions_discrete_mean=apply(predictions_discrete,c(1,2),mean)
      predictions_discrete_sd=apply(predictions_discrete,c(1,2),stats::sd)
      
      colnames(predictions_continuous_mean)=paste(alternatives,"continuous (mean)")
      colnames(predictions_continuous_sd)=paste(alternatives,"continuous (sd)")
      colnames(predictions_discrete_mean)=paste(alternatives,"discrete (mean)")
      colnames(predictions_discrete_sd)=paste(alternatives,"discrete (sd)")
      
      answer <- cbind(predictions_continuous_mean,predictions_continuous_sd,predictions_discrete_mean,predictions_discrete_sd)
      answer[!rows, ] <- NA
      return(answer)
    } else {
      cat("You have mixing")
      # check for intra-respondent, and if it exists, give a message, and average out third dimension. applies to all inputs (V,alpha,gamma,sigma)
      Ninter = apollo_inputs$apollo_draws$interNDraws

      ND=250
      consumption_overall=array(0,dim=c(nObs,length(alternatives),ND,Ninter))
      r=1
      while(r<=Ninter){
        Vdraw=V
        alphadraw=alpha
        gammadraw=gamma
        l=1
        while(l<=length(alternatives)){
          if(is.matrix(Vdraw[[l]])) Vdraw[[l]]=Vdraw[[l]][,r]
          if(is.matrix(alphadraw[[l]])) alphadraw[[l]]=alphadraw[[l]][,r]
          if(is.matrix(gammadraw[[l]])) gammadraw[[l]]=gammadraw[[l]][,r]
          l=l+1
        }
        consumption_overall[,,,r] = forecasting(Vdraw,alternatives,alphadraw,gammadraw,sigma,cost,avail,budget,continuousChoice)
        r=r+1
      }
      consumption <- apply(consumption_overall, MARGIN=c(1,2), mean) # Average intra-draws
      consumption[!rows, ] <- NA
      return( consumption )
    }
  }
}