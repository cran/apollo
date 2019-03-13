#' Calculates exploded logit probabilities
#'
#' Calculates the probabilities of an exploded logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' The function calculates the probability of a ranking as a product of logit models with gradually reducing availability, where scale differences can be allowed for.
#' @param el_settings List of inputs of the exploded logit model. It shoud contain the following.
#'                    \itemize{
#'                     \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                     \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                     \item choiceVars: List of numeric vectors. Contain choices for each position of the ranking. The list must be ordered with the best choice first, second best second, etc. It will usually be a list of columns from the database.
#'                     \item V: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                     \item scales: List of vectors. Scale factors of each logit model. Should have one element less than choiceVars. At least one element should be normalized to 1. If omitted, scale=1 for all positions is assumed.
#'                     \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                    }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate": Used for model estimation.
#'                        \item "prediction": Used for model predictions.
#'                        \item "validate": Used for validating input.
#'                        \item "zero_LL": Used for calculating null likelihood.
#'                        \item "conditionals": Used for calculating conditionals.
#'                        \item "output": Used for preparing output after model estimation.
#'                        \item "raw": Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": Not applicable.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @export
apollo_el <- function(el_settings, functionality){
  if(is.null(el_settings[["alternatives"]])) stop("The el_settings list needs to include an object called \"alternatives\"!")
  if(is.null(el_settings[["avail"]])) stop("The el_settings list needs to include an object called \"avail\"!")
  if(is.null(el_settings[["choiceVars"]])) stop("The el_settings list needs to include an object called \"choiceVars\"!")
  if(is.null(el_settings[["V"]])) stop("The el_settings list needs to include an object called \"V\"!")
  if(is.null(el_settings[["scales"]])) {
   el_settings[["scales"]]=vector(mode="list",length=length(el_settings[["alternatives"]])-1)
   el_settings[["scales"]] <- lapply(el_settings[["scales"]], function(a) 1)
  }
  if(is.null(el_settings[["rows"]])) el_settings[["rows"]]="all"

  alternatives = el_settings[["alternatives"]]
  avail        = el_settings[["avail"]]
  choiceVars   = el_settings[["choiceVars"]]
  V            = el_settings[["V"]]
  scales       = el_settings[["scales"]]
  rows         = el_settings[["rows"]]
  
  stages=length(choiceVars)
  
  if(functionality=="validate"){
    apollo_control <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_control,
                            error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    j=1
    while(j<=length(choiceVars)){
     choiceVars[[j]][!rows]=alternatives[1]
     j=j+1
    }

    if(length(rows)!=length(choiceVars[[1]])) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of the rows in the data!")
    if(length(scales)!=stages) stop("The object \"scales\" needs to either be the same length as the number of stages, i.e. one less than the number of alternatives!")

    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
      warning("Full availability of alternatives assumed for exploded logit component.")
    }
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]

    if(apollo_control$noValidation==FALSE){
      if(nAlts<3) stop("EL requires at least three alternatives")

      if(nObs==0) stop("No choices to model")

      choiceLabs <- unique(choiceVars[[1]])
      if(!all(altnames %in% names(V))) stop("Alternative labels in \"altnames\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternative labels in \"altnames\" do not match those in \"avail\".")

      if(!all(choiceLabs %in% altcodes)) stop("Value in choice column that is not included in altcodes.")

      avail_new=list()
      avail_new[[1]]=avail
      s=2
      while(s<=stages){
       avail_new[[s]]=avail_new[[s-1]]
       j=1
       while(j<=nAlts)
       {
        avail_new[[s]][j]=lapply(avail_new[[s]][j],"*",(choiceVars[[s-1]]!=j))
        j=j+1
       }
       s=s+1
      }
      
      s=1
      while(s<=stages){
        avail=avail_new[[s]]
        choiceVar=choiceVars[[s]]
        chosenunavail=0
        j=1
        while(j <= length(altnames)){
          if(sum((choiceVar==altcodes[j])*(avail[[j]]==0)*rows)) chosenunavail=1
          j=j+1
        }
        if(chosenunavail==1) stop("Some alternative(s) chosen despite being listed as unavailable in stage ",s,"\n")

        for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")
      s=s+1
      }
      
      cat("\nAll checks passed for exploded logit model component\n")

    }

    if(apollo_control$noDiagnostics==FALSE){
      if(avail_set==TRUE) warning("Availability not provided to 'apollo_el' (or some elements are NA). Full availability assumed.")

      s=1
      while(s<=stages){
      avail=avail_new[[s]]
      choiceVar=choiceVars[[s]]
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)

      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) 
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j] = sum(choiceVar==altcodes[j] & rows) 
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = altnames
      cat(paste("Overview of choices for exploded logit model component, stage ",s,":\n",sep=""))
      print(round(choicematrix,2))
      noAltAvail <- Reduce("+", avail_new[[s]])==0
      cat("There are ", sum(noAltAvail), " observation(s) with no available alternatives in stage ", s, ".\n", sep="")
      cat("\n")
      if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      s=s+1}
    }
  
    return(invisible(TRUE))
  }

 if(functionality=="zero_LL"){
    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    j=1
    while(j<=length(choiceVars)){
      choiceVars[[j]][!rows]=alternatives[1]
      j=j+1
    }
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }

    if(!anyNA(avail)) if(any(altnames != names(avail))) avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) 
    avail_new=list()
    avail_new[[1]]=avail
    nAvAlt <- rowSums(matrix(unlist(avail), ncol = length(avail))) 
    P = 1/nAvAlt 
    s=2
    while(s<=stages){
      avail_new[[s]]=avail_new[[s-1]]
      j=1
      while(j<=nAlts)
      {
        avail_new[[s]][j]=lapply(avail_new[[s]][j],"*",(choiceVars[[s-1]]!=j))
        j=j+1
      }
      nAvAlt <- rowSums(matrix(unlist(avail_new[[s]]), ncol = length(avail_new[[s]]))) 
      P = P*1/nAvAlt 
      s=s+1
    }
    P[!rows] = 1
    return(P)
  }

  if(functionality %in% c("estimate","conditionals","raw")){

    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    
    avail_new = list(avail)
    s=2
    while(s<=stages){
      avail_new[[s]] = mapply(function(a,j) a*(choiceVars[[s-1]]!=j), avail_new[[s-1]], as.list(1:nAlts), SIMPLIFY=FALSE)
      s=s+1
    }
    
    mnlSettings=list(
      alternatives = alternatives,
      avail        = avail_new[[1]],
      choiceVar    = choiceVars[[1]],
      V            = lapply(V,"*",scales[[1]]),
      rows         = rows
    )
    
    P = apollo_mnl(mnlSettings, functionality) 
    s=2
    while(s<=stages){
      mnlSettings$avail     = avail_new[[s]]
      mnlSettings$choiceVar = choiceVars[[s]]
      mnlSettings$V         = lapply(V,"*",scales[[s]])
      mnlSettings$rows      = mnlSettings$rows & (Reduce("+", avail_new[[s]])>0)
      P=P*apollo_mnl(mnlSettings, functionality)
      s=s+1
    }
    
    return(P)
  }

  if(functionality=="prediction"){
   return(NA)
  }
    
  if(functionality=="output"){

    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    
    P = apollo_el(el_settings, functionality="estimate")
    
    
    avail_new = list(avail)
    for(s in 2:stages){
      avail_new[[s]] = avail_new[[s-1]]
      for(j in 1:nAlts) avail_new[[s]][j] = lapply(avail_new[[s]][j],"*",(choiceVars[[s-1]]!=j))
    }
    
    s=1
    while(s<=stages){
      avail     = avail_new[[s]]
      choiceVar = choiceVars[[s]]
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      availprint   = colSums(rows*matrix(unlist(avail), ncol = length(avail))) 
      choicematrix = matrix(0,nrow=4,ncol=nAlts)
      choicematrix[1,] = availprint
      for(j in 1:nAlts) choicematrix[2,j] = sum(choiceVar==alternatives[j] & rows) 
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = names(alternatives)
      
      apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
                                  error=function(e){
                                    cat("apollo_el could not retrieve apollo_control. No diagnostics in output.\n")
                                    return(NA)
                                  } )
      if(!(length(apollo_control)==1 && is.na(apollo_control))){
        fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
        fileName <- file.path(tempdir(),fileName)
        fileConn <- tryCatch( file(fileName, open="at"),
                              error=function(e){
                                cat('apollo_el could not write diagnostics to temporary file. No diagnostics in output.\n')
                                return(NA)
                              })
        if(!anyNA(fileConn)){
          sink(fileConn)
          on.exit({if(sink.number()>0) sink(); close(fileConn)})
          if(apollo_control$noDiagnostics==FALSE){
            cat(paste("Overview of choices for exploded logit model component, stage ",s,":\n",sep=""))
            print(round(choicematrix,2))
            noAltAvail <- Reduce("+", avail_new[[s]])==0
            cat("There are ", sum(noAltAvail), " observation(s) with no available alternatives in stage ", s, ".\n", sep="")
            cat("\n")}
          if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
          if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
        }
      }      
  
      s=s+1}
 
     return(P)
  }
}

