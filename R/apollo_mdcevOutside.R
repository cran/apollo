#' Calculates MDCEV likelihoods with an outside good.
#' 
#' Calculates the likelihood of a Multiple Discrete Continuous Extreme Value (MDCEV) model with an outside good.
#' 
#' @param V Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#' @param alternatives Character vector. Names of alternatives, elements must match the names in list 'V'.
#' @param alpha Named list. Alpha parameters for each alternative, including for the outside good. As many elements as alternatives.
#' @param gamma Named list. Gamma parameters for each alternative, excluding the outside good. As many elements as inside good alternatives.
#' @param sigma Numeric scalar. Scale parameter of the model extreme value type I error.
#' @param cost Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#' @param avail Named list. Availabilities of alternatives, one element per alternative. Names of elements must match those in argument 'alternatives'. Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation. If all alternatives are always available, then user can just omit this argument.
#' @param continuousChoice Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#' @param budget  Numeric vector. Budget for each observation.
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
#' @param minConsumption  Named list of scalars or numeric vectors. Minimum consumption of the alternatives, if consumed. As many elements as alternatives. Names must match those in \code{alternatives}.
#' @param outsideName Character. Alternative name for the outside good. The first good is assumed to be the outside one. Default is "outside"
#' @param rows  Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#' @param componentName Character. Name given to model component.
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
apollo_mdcevOutside <- function(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice, budget, functionality, minConsumption=NA, outsideName="outside", rows="all", componentName="MDCEV"){
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # This functionality should be dealed with by apollo_mdcev.R
    return(TRUE)
  }

  # ############################# #
  #### functionality="zero_LL" ####
  # ############################# #
  
  if(functionality=="zero_LL"){
    nObs  <- length(continuousChoice[[1]])
    P <- rep(NA,nObs)
    P[!rows] <- 1
    return(P)
  }
  
  # ############################################### #
  #### functionality="estimate/conditionals/raw" ####
  # ################################################# #
  
  if(functionality %in% c("estimate","conditionals","raw")){
    
    # Store useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    discrete_choice <- lapply(continuousChoice, function(x) x>0)
    
    # Compute V
    V[[1]]=(alpha[[1]]-1)*log(continuousChoice[[1]])
    j=2
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
    
    # First term
    term1=(1-totalChosen)*log(sigma)
    
    # Second term
    logfi=list()
    j=2
    while(j<=length(V)){
      if(!anyNA(minConsumption)){
        tmp <- continuousChoice[[j]]-(continuousChoice[[j]]>=minConsumption[[j]])*minConsumption[[j]]
        logfi[[j-1]]=avail[[j]]*( log(1-alpha[[j]]) - log(tmp + gamma[[j]]) )
      } else {
        logfi[[j-1]]=avail[[j]]*(log(1-alpha[[j]])-log(continuousChoice[[j]]+gamma[[j]]))
      }
      j=j+1
    }
    
    term2=log(1-alpha[[1]])-log(continuousChoice[[1]])
    j=2
    while(j<=length(V)){
      term2=term2+avail[[j]]*(logfi[[j-1]]*discrete_choice[[j]])
      j=j+1
    }
    
    # Third term
    term3=continuousChoice[[1]]/(1-alpha[[1]])
    j=2
    while(j<=length(V)){
      term3=term3+avail[[j]]*(cost[[j]]/exp(logfi[[j-1]]) * discrete_choice[[j]])
      j=j+1}
    term3 = log(term3)
    
    # Fourth term
    term4_1 = V[[1]]/sigma
    term4_2 = exp(V[[1]]/sigma)
    j=2
    while(j<=length(V)){
      term4_1 = term4_1+avail[[j]]*(V[[j]]/sigma * discrete_choice[[j]])
      term4_2 = term4_2+avail[[j]]*exp(V[[j]]/sigma)
      j=j+1
    }
    term4_2 = totalChosen * log(term4_2) 
    term4 =  term4_1 - term4_2 
    
    # Fifth term: log of factorial
    term5 = lfactorial(totalChosen-1)
    
    # probability is simply the exp of the sum of the logs of the individual terms
    P = exp(term1 + term2 + term3 + term4 + term5)
    if(is.vector(P)) P[!rows]   <- 1
    if(is.matrix(P)) P[!rows,]  <- 1
    if(is.array(P) && length(dim(P))==3) P[!rows,,] <- 1
    
    # make the chosen unavailable alternatives have a likelihood of zero
    choseUnavail <- mapply(function(m,a) m & !a, discrete_choice, avail, SIMPLIFY=TRUE)
    choseUnavail <- rowSums(choseUnavail)>0
    if(is.vector(P)) P[choseUnavail]   <- 0
    if(is.matrix(P)) P[choseUnavail,]  <- 0
    if(is.array(P) && length(dim(P))==3) P[choseUnavail,,] <- 0
    
    return(P)
    
  }
  
  # ############################## #
  #### functionality="output" ####
  # ############################## #
  
  if(functionality=="output"){
    # Calculate likelihood
    P <- apollo_mdcevOutside(V, alternatives, alpha, gamma, sigma, cost, avail, continuousChoice, budget, functionality="estimate", minConsumption, rows)
    
    # Useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    discrete_choice <- lapply(continuousChoice, function(x) x>0)
    
    # turn scalar availabilities into vectors
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
    
    # Calculate data summary
    choicematrix = matrix(0,nrow=4,ncol=length(V))
    choicematrix[1,] = colSums(matrix(unlist(avail), ncol = length(avail))[rows,])
    for(j in 1:length(V)){
      choicematrix[2,j]=sum(discrete_choice[[j]][rows]==1)
      choicematrix[3,j]=sum(continuousChoice[[j]][rows])/choicematrix[1,j]
      choicematrix[4,j]=sum(continuousChoice[[j]][rows])/choicematrix[2,j]
    }
    choicematrix <- apply(choicematrix, MARGIN=2, function(x) {x[!is.finite(x)] <- 0; return(x)})
    rownames(choicematrix) = c("Times available","Observations in which chosen","Average consumption when available","Average consumption when chosen")
    colnames(choicematrix) = names(V)
    
    # Write to apolloLog
    content <- list(round(choicematrix,2))
    if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog(title   = paste0("Overview of choices  for model component \"",componentName,"\":"), 
                  content = content, apolloLog)
    apollo_reportModelTypeLog(modelType="MDCEV", apolloLog)
    
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #
  
  if(functionality=="prediction"){
    
    # Store useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    
    # first include a check to make sure that all alphas are equal across alts
    equality_check=1
    j=2
    while(j<= length(alpha)){
      if(alpha[[j]]!=alpha[[1]]) equality_check=0
      j=j+1
    }
    if(equality_check!=1) stop("\nMDCEV prediction for model component \"",componentName,"\" only implemented for profile where alpha is generic across all products!")
    
    # include a check to make sure there are no minimum consumption limits
    if(!anyNA(minConsumption)) stop("\nMDCEV prediction for model component \"",componentName,"\" only implemented for case without non-zero minimum consumptions!")
    
    forecasting=function(V,alternatives,alpha,gamma,sigma,cost,avail,budget,continuousChoice)
    {
      set.seed(99)
      
      ND=250
      nObs=length(continuousChoice[[1]])[1]
      consumption_total=array(0,dim=c(nObs,length(alternatives),ND))
      
      # create uniform draws for value draws
      # this needs to change still if we also have random alpha, gamma or V
      
      Ndraws=apollo_mlhs(ND,length(alternatives),nObs)
      draws_base=cbind(x=rep(seq(1,ND),nObs),Ndraws)
      
      k=1
      cat("\n0%")
      while(k<(ND+1))
      {
        draws=subset(draws_base,draws_base[,1]==k)
        draws=draws[,2:ncol(draws)]
        
        draws=sigma*(-log(-log(draws)))
        
        V_trans=V
        
        V_trans[[1]]=exp(V_trans[[1]]+draws[,1])
        
        j=2
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
        
        # need a check here to see if we're using a given column for example of the cube
        
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
        orderofV=rank(-V_trans[i,2:NALT])
        M=1
        stopping=0
        use=(orderofV<M)
        while(stopping<1)
        {
          #step2
          lambda_1=((budget[i]+sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use)))
          lambda_21=(V_trans[i,1])^(1/(1-alpha[i,1]))
          lambda_22=sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use*(V_trans[i,2:NALT])^(1/(1-alpha[i,1])))
          lambda=(lambda_1/(lambda_21+lambda_22))^(alpha[i,1]-1)
          if((sum(lambda<V_trans[i,2:NALT]))<(M))
          {              
            #step3
            consumption_outside_1=(V_trans[i,1])^(1/(1-alpha[i,1]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
            consumption_outside_2=lambda_21+lambda_22
            consumption[i,1]=consumption_outside_1/consumption_outside_2
            consumption_inside_1=(V_trans[i,2:NALT])^(1/(1-alpha[i,2:NALT]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
            consumption_inside_2=consumption_outside_2
            consumption[i,2:NALT]=use*((consumption_inside_1/consumption_inside_2)-1)*gamma[i,(2:ncol(gamma))]
            stopping=1
          }
          else
          {
            #step4
            M=M+1  
            use=(orderofV<M)
            if(M==(sum(avail[i,2:NALT])+1))
            {
              lambda_22=sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use*(V_trans[i,2:NALT])^(1/(1-alpha[i,1])))
              consumption_outside_1=(V_trans[i,1])^(1/(1-alpha[i,1]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
              consumption_outside_2=lambda_21+lambda_22
              consumption[i,1]=consumption_outside_1/consumption_outside_2
              consumption_inside_1=(V_trans[i,2:NALT])^(1/(1-alpha[i,2:NALT]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
              consumption_inside_2=consumption_outside_2
              consumption[i,2:NALT]=use*(consumption_inside_1/consumption_inside_2-1)*gamma[i,(2:ncol(gamma))]
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
      
      cat("\nNow producing forecasts for model component \"",componentName,"\".",sep="")
      cat("\nA matrix with one row per observation and the following columns will be returned:")
      cat("\n1. Means of predicted continuous consumptions (one column per product)")
      cat("\n2. Std err of predicted continuous consumptions (one column per product)")
      cat("\n3. Means of predicted discrete consumptions (one column per product)")
      cat("\n4. Std err of predicted discrete consumptions (one column per product)")
      cat("\n\nThis may take a while!")
      
      consumption=forecasting(V,alternatives,alpha,gamma,sigma,cost,avail,budget,continuousChoice)
      
      dimnames(consumption)[[2]]=c(outsideName, alternatives[2:length(alternatives)])
      
      predictions_continuous_mean = apply(consumption,c(1,2),mean)
      predictions_continuous_sd   = apply(consumption,c(1,2),stats::sd)
      predictions_discrete        = (consumption>0)
      predictions_discrete_mean   = apply(predictions_discrete,c(1,2),mean)
      predictions_discrete_sd     = apply(predictions_discrete,c(1,2),stats::sd)
      
      tmp <- c(outsideName, alternatives[2:length(alternatives)])
      colnames(predictions_continuous_mean) = paste(tmp,"continuous (mean)")
      colnames(predictions_continuous_sd)   = paste(tmp,"continuous (sd)")
      colnames(predictions_discrete_mean)   = paste(tmp,"discrete (mean)")
      colnames(predictions_discrete_sd)     = paste(tmp,"discrete (sd)")
      
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
