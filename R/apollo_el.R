#' Calculates exploded logit probabilities
#'
#' Calculates the probabilities of an exploded logit model and can also perform other operations based on the value of the \code{functionality} argument.
#' The function calculates the probability of a ranking as a product of logit models with gradually reducing availability, where scale differences can be allowed for.
#' @param el_settings List of inputs of the exploded logit model. It shoud contain the following.
#'                    \itemize{
#'                     \item \strong{\code{"alternatives"}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                     \item \strong{\code{"avail"}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                     \item \strong{\code{"choiceVars"}}: List of numeric vectors. Contain choices for each position of the ranking. The list must be ordered with the best choice first, second best second, etc. It will usually be a list of columns from the database.
#'                     \item \strong{\code{"V"}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                     \item \strong{\code{"scales"}}: List of vectors. Scale factors of each logit model. Should have one element less than choiceVars. At least one element should be normalized to 1. If omitted, scale=1 for all positions is assumed.
#'                     \item \strong{\code{"rows"}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                     \item \strong{\code{"componentName"}}: Character. Name given to model component.
#'                    }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"}: Used for model estimation.
#'                        \item \code{"prediction"}: Used for model predictions.
#'                        \item \code{"validate"}: Used for validating input.
#'                        \item \code{"zero_LL"}: Used for calculating null likelihood.
#'                        \item \code{"conditionals"}: Used for calculating conditionals.
#'                        \item \code{"output"}: Used for preparing output after model estimation.
#'                        \item \code{"raw"}: Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: Not applicable (\code{NA}).
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'         }
#' @importFrom stats setNames
#' @export
apollo_el <- function(el_settings, functionality){
  if(is.null(el_settings[["componentName"]])       ) el_settings[["componentName"]]="EL"
  componentName= el_settings[["componentName"]]
  
  if(is.null(el_settings[["alternatives"]])) stop("The el_settings list for model component \"",componentName,"\" needs to include an object called \"alternatives\"!")
  if(is.null(el_settings[["avail"]])) stop("The el_settings list for model component \"",componentName,"\" needs to include an object called \"avail\"!")
  if(is.null(el_settings[["choiceVars"]])) stop("The el_settings list for model component \"",componentName,"\" needs to include an object called \"choiceVars\"!")
  if(is.null(el_settings[["V"]])) stop("The el_settings list for model component \"",componentName,"\" needs to include an object called \"V\"!")
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
  nObs <- tryCatch(nrow( get("apollo_inputs", parent.frame(), inherits=FALSE)$database ),
                   error=function(e){
                     lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenA <- sapply(avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenC <- sapply(choiceVars, length)
                     return(max(lenV, lenA, lenC))
                   })
  
  ### Format checks
  # alternatives
  test <- is.vector(alternatives) & !is.null(names(alternatives))
  if(!test) stop("The \"alternatives\" argument for model component \"",componentName,"\" needs to be a named vector")
  # avail
  test <- is.list(avail) || (length(avail)==1 && avail==1)
  if(!test) stop("The \"avail\" argument for model component \"",componentName,"\" needs to be a list or set to 1")
  if(is.list(avail)){
    lenA <- sapply(avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
    test <- all(lenA==nObs | lenA==1)
    if(!test) stop("All entries in \"avail\" for model component \"",componentName,"\" need to be a scalar or a vector with one entry per observation in the \"database\"")
  }
  # choiceVar
  test <- is.list(choiceVars) && (stages==length(V)-1)
  if(!test) stop("The \"choiceVars\" argument for model component \"",componentName,"\" needs to be a list, with one vector per stage, each with one entry per observation in the \"database\"")
  # V
  if(!is.list(V)) stop("The \"V\" argument for model component \"",componentName,"\" needs to be a list")
  lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
  test <- all(lenV==nObs | lenV==1)
  if(!test) stop("Each element of \"V\" for model component \"",componentName,"\" needs to be a scalar or a vector/matrix/cube with one row per observation in the \"database\"")  
  if(length(scales)!=stages) stop("The object \"scales\" for model component \"",componentName,"\" needs to be the same length as the number of stages, i.e. one less than the number of alternatives!")
  # rows
  test <- is.vector(rows) && ( (is.logical(rows) && length(rows)==nObs) || (length(rows)==1 && rows=="all") )
  if(!test) stop("The \"rows\" argument for model component \"",componentName,"\" needs to be \"all\" or a vector of boolean statements with one entry per observation in the \"database\"")
  # functionality
  test <- functionality %in% c("estimate","prediction","validate","zero_LL","conditionals","output","raw")
  if(!test) stop("Non-permissable setting for \"functionality\" for model component \"",componentName,"\"")
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Store useful values
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
    
    # check rows statement
    
    # Check that alternatives are named in altcodes and V
    if(is.null(altnames) || is.null(altcodes) || is.null(names(V))) stop("Alternatives for model component \"",componentName,"\" must be named, both in 'alternatives' and 'V'.")
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
      #warning("Full availability of alternatives assumed for exploded logit component.")
    }
    # Reorder V to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]
    
    if(apollo_control$noValidation==FALSE){
      # Check that there are at least two alternatives
      if(nAlts<3) stop("Model component \"",componentName,"\"  requires at least three alternatives")
      
      # Check that choice vector is not empty
      if(nObs==0) stop("No data for model component \"",componentName,"\"")
      
      # Check that labels in choice match those in the utilities and availabilities
      choiceLabs <- unique(choiceVars[[1]])
      if(!all(altnames %in% names(V))) stop("The names of the alternatives for model component \"",componentName,"\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("The names of the alternatives for model component \"",componentName,"\" do not match those in \"avail\".")
      
      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop("The data contains values for \"choiceVar\" for model component \"",componentName,"\" that are not included in \"alternatives\".")
      
      # create new availabilities list, with one list per stage
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
      
      # check that nothing unavailable is chosen
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
        if(chosenunavail==1) stop("Some alternative(s) for model component \"",componentName,"\" chosen despite being listed as unavailable in stage ",s,"\n")
        
        # check that all availabilities are either 0 or 1
        for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values for model component \"",componentName,"\" are not 0 or 1.")
        s=s+1
      }
      
      # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
      V <- mapply(function(v,a) apollo_setRows(v, !a | !rows, 0), V, avail, SIMPLIFY=FALSE)
      if(any(sapply(V, anyNA))) cat("\nAt least one utility for model component \"",componentName,"\" contains one or more NA values")
      
      
    }
    
    if(apollo_control$noDiagnostics==FALSE){
      #if(avail_set==TRUE) warning("Availability not provided to 'apollo_el' (or some elements are NA). Full availability assumed.")
      
      content <- list()
      s=1
      while(s<=stages){
        avail=avail_new[[s]]
        choiceVar=choiceVars[[s]]
        # turn scalar availabilities into vectors
        for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
        
        # Construct summary table of availabilities and market share
        availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
        choicematrix = matrix(0,nrow=4,ncol=length(altnames))
        choicematrix[1,] = availprint
        j=1
        while(j<= length(altnames)){
          choicematrix[2,j] = sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
          j=j+1
        }
        choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
        choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
        choicematrix[4,!is.finite(choicematrix[4,])] <- 0
        rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
        colnames(choicematrix) = altnames
        #cat(paste("Overview of choices for exploded logit model component, stage ",s,":\n",sep=""))
        #print(round(choicematrix,2))
        content[[length(content) + 1]] <- paste("Overview of choices for model component \"",componentName,"\", stage ",s,":\n",sep="")
        content[[length(content) + 1]] <- round(choicematrix,2)
        
        # Check if there are shorter rankings
        noAltAvail <- sum(Reduce("+", avail_new[[s]])==0)
        if (noAltAvail>0) cat("There are ", noAltAvail, " observation(s) with no available alternatives in stage ", s, ".\n", sep="")
        #cat("\n")
        if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
        if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
        s=s+1
      }
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog(title   = paste0("Overview of choices  for model component \"",componentName,"\":"), 
                    content = content, apolloLog)
    }
    
    testL=apollo_el(el_settings, functionality="estimate")
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"",componentName,"\"")
    if(any(testL==0)) cat("\nSome observations have zero probability at starting value for model component \"",componentName,"\"")
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # Store useful values
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
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }
    
    if(!anyNA(avail)) if(any(altnames != names(avail))) avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    # create new availabilities list, with one list per stage
    avail_new=list()
    avail_new[[1]]=avail
    nAvAlt <- rowSums(matrix(unlist(avail), ncol = length(avail))) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    s=2
    while(s<=stages){
      avail_new[[s]]=avail_new[[s-1]]
      j=1
      while(j<=nAlts)
      {
        avail_new[[s]][j]=lapply(avail_new[[s]][j],"*",(choiceVars[[s-1]]!=j))
        j=j+1
      }
      nAvAlt <- rowSums(matrix(unlist(avail_new[[s]]), ncol = length(avail_new[[s]]))) # number of available alts in each observation
      P = P*1/nAvAlt # likelihood at zero
      s=s+1
    }
    P[!rows] = 1
    return(P)
  }
  
  # ############################################################ #
  #### functionality="estimate/conditionals/raw"              ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","conditionals","raw")){
    
    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    avail_set <- FALSE
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), names(alternatives))
    }
    
    # create new availabilities list, with one list per stage
    avail_new = list(avail)
    s=2
    while(s<=stages){
      avail_new[[s]] = mapply(function(a,j) a*(choiceVars[[s-1]]!=j), avail_new[[s-1]], as.list(1:nAlts), SIMPLIFY=FALSE)
      s=s+1
    }
    
    mnl_settings=list(
       alternatives = alternatives,
       avail        = avail_new[[1]],
       choiceVar    = choiceVars[[1]],
       V            = lapply(V,"*",scales[[1]]),
       rows         = rows
     )
    
    # Calculate MNL for each stage
    P = apollo_mnl(mnl_settings, functionality) # 1st stage
    s=2
    while(s<=stages){
      # Update settings
      mnl_settings$avail     = avail_new[[s]]
      mnl_settings$choiceVar = choiceVars[[s]]
      mnl_settings$V         = lapply(V,"*",scales[[s]])
      mnl_settings$rows      = mnl_settings$rows & (Reduce("+", avail_new[[s]])>0)
      # Calculate probabilities
      P=P*apollo_mnl(mnl_settings, functionality)
      s=s+1
    }
    
    return(P)
  }
  
  # ################################ #
  #### functionality="prediction" ####
  # ################################ #
  
  if(functionality=="prediction"){
    return(NA)
  }
  
  # ############################## #
  #### functionality="output" ####
  # ############################## #
  
  if(functionality=="output"){
    
    # Store useful values
    nObs  <- length(choiceVars[[1]])
    nAlts <- length(V)
    altnames <- names(alternatives)
    altcodes <- alternatives
    avail_set <- FALSE
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVars[[1]]))
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), names(alternatives))
    }
    
    P = apollo_el(el_settings, functionality="estimate")
    
    
    # create new availabilities list, with one list per stage
    avail_new = list(avail)
    for(s in 2:stages){
      avail_new[[s]] = avail_new[[s-1]]
      for(j in 1:nAlts) avail_new[[s]][j] = lapply(avail_new[[s]][j],"*",(choiceVars[[s-1]]!=j))
    }
    
    ### Write to apolloLog
    content <- list()
    s=1
    while(s<=stages){
      avail=avail_new[[s]]
      choiceVar=choiceVars[[s]]
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j] = sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = altnames
      #cat(paste("Overview of choices for exploded logit model component, stage ",s,":\n",sep=""))
      #print(round(choicematrix,2))
      content[[length(content) + 1]] <- paste("Overview of choices for model component \"",componentName,"\", stage ",s,":\n",sep="")
      content[[length(content) + 1]] <- round(choicematrix,2)
      
      # Check if there are shorter rankings
      noAltAvail <- sum(Reduce("+", avail_new[[s]])==0)
      if (noAltAvail>0) cat("There are ", noAltAvail, " observation(s) with no available alternatives in stage ", s, ".\n", sep="")
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!" #cat("Warning: some alternatives are never chosen in your data!\n")
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!" #cat("Warning: some alternatives are always chosen when available!\n")
      s=s+1
    }
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog(title   = paste0("Overview of choices  for model component \"",componentName,"\":"), 
                  content = content, apolloLog)
    apollo_reportModelTypeLog(modelType="EL", apolloLog)
    
    return(P)
  }
}

