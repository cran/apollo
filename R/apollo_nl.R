#' Calculates probabilities of a nested logit
#'
#' Calculates probabilities of a nested logit model.
#'
#' In this implementation of the nested logit model, each nest must have a lambda parameter associated to it.
#' For the model to be consistent with utility maximisation, the estimated value of the Lambda parameter of all nests
#' should be between 0 and 1. Lambda parameters are inversely proportional to the correlation between the error terms of 
#' alternatives in a nest. If lambda=1, then there is no relevant correlation between the unobserved
#' utility of alternatives in that nest.
#' The tree must contain an upper nest called \code{"root"}. The lambda parameter of the root is automatically
#' set to 1 if not specified in \code{nlNests}. And while setting it to another value is possible, it is not
#' recommended.
#' @param nl_settings List of inputs of the NL model. It shoud contain the following.
#'                    \itemize{
#'                       \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item choiceVar: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item V: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item nlNests: List of numeric scalars or vectors. Lambda parameters for each nest. Elements must be named with the nest name. The lambda at the root is fixed to 1 if excluded (recommended).
#'                       \item nlStructure: Named list of character vectors. As many elements as nests, it must include the "root". Each element contains the names of the nests or alternatives that belong to it. Element names must match those in \code{nlNests}.
#'                       \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
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
#'           \item "prediction": List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @examples
#' ### Load data
#' data(apollo_modeChoiceData)
#' database <- apollo_modeChoiceData
#' rm(apollo_modeChoiceData)
#'
#' ### Parameters
#' b = list(asc_1=0, asc_2=0, asc_3=0, asc_4=0, tt=0, tc=0, acc=0, lambda=0.5)
#'
#' V = list()
#' V[['car' ]] = b$asc_1 + b$tt*database$time_car  + b$tc*database$cost_car
#' V[['bus' ]] = b$asc_2 + b$tt*database$time_bus  + b$tc*database$cost_bus  +
#'               b$acc*database$access_bus
#' V[['air' ]] = b$asc_3 + b$tt*database$time_air  + b$tc*database$cost_air  +
#'               b$acc*database$access_air
#' V[['rail']] = b$asc_4 + b$tt*database$time_rail + b$tc*database$cost_rail +
#'               b$acc*database$access_rail
#'
#' ### NL settings
#' nl_settings <- list(
#'    alternatives = c(car=1, bus=2, air=3, rail=4),
#'    avail        = list(car=database$av_car, bus=database$av_bus,
#'                        air=database$av_air, rail=database$av_rail),
#'    choiceVar    = database$choice,
#'    V            = V,
#'    nlNests      = list(root=1, public=b$lambda),
#'    nlStructure  = list(root=c("car", "public"), public=c("bus","air","rail"))
#' )
#'
#' ### Compute choice probabilities using NL model
#' apollo_nl(nl_settings, functionality="estimate")
#' @export
apollo_nl <- function(nl_settings, functionality){
  if(is.null(nl_settings[["alternatives"]])) stop("The nl_settings list needs to include an object called \"alternatives\"!")
  if(is.null(nl_settings[["avail"]])) stop("The nl_settings list needs to include an object called \"avail\"!")
  if(is.null(nl_settings[["choiceVar"]])) stop("The nl_settings list needs to include an object called \"choiceVar\"!")
  if(is.null(nl_settings[["V"]])) stop("The nl_settings list needs to include an object called \"V\"!")
  if(is.null(nl_settings[["nlNests"]])) stop("The nl_settings list needs to include an object called \"nlNests\"!")
  if(is.null(nl_settings[["nlStructure"]])) stop("The nl_settings list needs to include an object called \"nlStructure\"!")
  if(is.null(nl_settings[["rows"]])) nl_settings[["rows"]]="all"
  
  alternatives = nl_settings[["alternatives"]]
  avail        = nl_settings[["avail"]]
  choiceVar    = nl_settings[["choiceVar"]]
  V            = nl_settings[["V"]]
  rows         = nl_settings[["rows"]]
  nlNests      = nl_settings[["nlNests"]]
  nlStructure  = nl_settings[["nlStructure"]]
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Store useful values
    apollo_control <- tryCatch(get("apollo_control", parent.frame(), inherits=FALSE),
                               error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    root_set  <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    nestnames <- names(nlStructure)
    if(!("root" %in% names(nlNests))){
      root_set <- TRUE
      nlNests["root"] <- 1
      #warning("Root lambda parameter set to 1.")
    }
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # check rows statement
    if(length(rows)!=length(choiceVar)) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of the rows in the data!")
    
    # Check that alternatives are named in altcodes and V
    if(is.null(altnames) || is.null(altcodes) || is.null(names(V))) stop("Alternatives must be named, both in 'altnernatives' and 'V'.")
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
      #warning("Full availability of alternatives assumed for NL component.")
    }
    # Reorder V to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(any(altnames != names(avail))) avail <- avail[altnames]
    
    if(apollo_control$noValidation==FALSE){
      
      # Check there are at least three alternatives
      if(nAlts<3) stop("NL requires at least three alternatives")
      
      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")
      
      # Check that labels in choice match those in the utilities and availabilities
      choiceLabs <- unique(choiceVar)
      if(!all(altnames %in% names(V))) stop("Alternative labels in \"altnames\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternative labels in \"altnames\" do not match those in \"avail\".")
      
      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop("Value in choice column that is not included in altcodes.")
      
      # check that nothing unavailable is chosen
      chosenunavail=0
      j=1
      while(j <= length(altnames)){
        if(sum((choiceVar==altcodes[j])*(avail[[j]]==0)*rows)) chosenunavail=1
        j=j+1
      }
      if(chosenunavail==1) stop("Some alternative(s) chosen despite being listed as unavailable\n")
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")
      
      
      # checks that are specific to nlStructure component
      
      allElements <- c("root", unlist(nlStructure))
      if(nlNests["root"]!=1) stop("The root lambda parameter should be equal to 1.")
      if( !all(altnames %in% allElements) ) stop("All alternatives must be included in the tree structure.")
      if( !all(nestnames %in% allElements) ) stop("All nests must be included in the tree structure.")
      if( (length(nestnames)+length(altnames))!=length(allElements) ) stop("Tree structure is inconsistent. Each element must appear only once.")
      if( !all(names(nlNests) %in% names(nlStructure)) | !all(names(nlStructure) %in% names(nlNests)) ) stop("All nests in argument 'nlNests' should be in 'nlStructure', and vice versa (including 'root').")
      
      if(is.null(nlStructure[["root"]])) stop("Tree structure is missing an element called root!")
      combined_elements="root"
      j=1
      while(j<= length(nlStructure)){
        combined_elements=c(combined_elements,nlStructure[[j]])
        j=j+1
      }
      
      j=1
      while(j<= length(altnames)){
        if(sum(nestnames==altnames[j])) stop("A nest cannot have the same name as an alternative!")
        if(sum(combined_elements==altnames[j])!=1) stop("An alternative needs to appear exactly once in a tree!")
        j=j+1
      }
      
      j=1
      while(j<= length(nlStructure)){
        if(sum(nestnames==names(nlStructure)[j])!=1) stop("A defined nest needs to appear exactly once in a tree!")
        j=j+1
      }
      
      j=1
      while(j<= length(nestnames)){
        if(sum(altnames==nestnames[j])) stop("A nest cannot have the same name as an alternative!")
        if(sum(combined_elements==nestnames[j])!=1) stop("A defined nest needs to appear exactly once in a tree!")
        j=j+1
      }
      
      # reorder nlStructure if needed
      nlStructure_ordered=list()
      element_list="root"
      j=1
      while(j>0){
        k=1
        temp=rep(TRUE,length(element_list))
        while(k<= length(element_list)){
          if(element_list[k] %in% altnames) temp[k]=FALSE
          k=k+1
        }
        element_list=element_list[temp]
        if(length(element_list)>0){
          nlStructure_ordered[[element_list[1]]]=nlStructure[[element_list[1]]]
          element_list=c(element_list,nlStructure[[element_list[1]]])
          element_list=element_list[-1]
        }
        j=length(element_list)
      }
      
      nlStructure=nlStructure_ordered
      
      # confirm all checks are passed, print output and return TRUE
      #cat("\nAll checks passed for NL model component\n")
    }
    if(apollo_control$noDiagnostics==FALSE){
      
      #if(avail_set==TRUE) warning("Availability not provided to 'apollo_nl' (or some elements are NA). Full availability assumed.")
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
      colnames(choicematrix) = altnames
      
      # Order tree structure
      nlStructure_ordered=list()
      element_list="root"
      j=1
      while(j>0){
        k=1
        temp=rep(TRUE,length(element_list))
        while(k<= length(element_list)){
          if(element_list[k] %in% altnames) temp[k]=FALSE
          k=k+1
        }
        element_list=element_list[temp]
        if(length(element_list)>0){
          nlStructure_ordered[[element_list[1]]]=nlStructure[[element_list[1]]]
          element_list=c(element_list,nlStructure[[element_list[1]]])
          element_list=element_list[-1]
        }
        j=length(element_list)
      }
      nlStructure=nlStructure_ordered
      
      # Calculate ancestors
      ancestors=list()
      j=1
      while(j<= length(altnames)){
        altJ <- altnames[[j]]
        ancestors[[altJ]] = altJ
        current = altJ
        k = length(nlStructure)
        while(k>0){
          if(current %in% nlStructure[[k]]){
            ancestors[[altnames[[j]]]] = c(ancestors[[altJ]],names(nlStructure)[k])
            current = names(nlStructure)[k]
          }
          k=k-1
        }
        j=j+1
      }
      
      # Load function to print tree structure and print tree structure
      print_tree=function(nlStructure, ancestors){
        
        print_tree_level = function(nlStructure, component, preceding_nest_layer, space){
          j=1
          if(preceding_nest_layer!=0) space=c(space,"  |")
          while(j<=length(nlStructure[[component]])){
            space <- gsub("[']", " ", space)
            if(j==length(nlStructure[[component]])) space[length(space)] <- gsub("[|]", "'", space[length(space)])
            if(nlStructure[[component]][j] %in% altnames){
              depth <- length(space)
              cat("\n",space,rep("-",3*(maxDepth-depth)),"-Alternative: ",nlStructure[[component]][j], sep="")
            } else {
              #cat("\n", space, "-Nest: ", nlStructure[[component]][j], sep="")
              cat("\n",space,"-Nest: ",nlStructure[[component]][j]," (",round(nlNests[[nlStructure[[component]][j]]],4), ")", sep="")
              print_tree_level(nlStructure, nlStructure[[component]][j], preceding_nest_layer+1, space)
            }
            j=j+1
          }
        }
        
        maxDepth <- max(sapply(ancestors, length))-1
        #cat("\nNest: ",names(nlStructure)[[1]], sep="")
        cat("Nest: ",names(nlStructure)[[1]]," (",round(nlNests[[names(nlStructure)[[1]]]],4),")", sep="")
        
        print_tree_level(nlStructure, "root", preceding_nest_layer=0, space="|")
      }
      
      
      ### Printing diagnostics
      #cat('Overview of choices for NL model component:\n')
      #print(round(choicematrix,2))
      #cat("\n")
      #if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      #if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      #cat("Nested logit structure:")
      #print_tree(nlStructure, ancestors)
      #cat("\n")
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content)+1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content)+1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      if(root_set==TRUE) content[[length(content)+1]] <- "Notice: Root lambda parameter set to 1."
      content[[length(content)+1]] <- "Nested structure:"
      content[[length(content)+1]] <- capture.output(print_tree(nlStructure, ancestors))
      #content <- list(round(choicematrix,2),
      #                ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
      #                ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""),
      #                "Nested structure:",
      #                capture.output(print_tree(nlStructure, ancestors)))
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for NL model component:", content, apolloLog)
    }
    
    return(TRUE)
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    nestnames <- names(nlStructure)
    if(!("root" %in% names(nlNests))) nlNests["root"] <- 1
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }
    
    if(!anyNA(avail)) if(any(altnames != names(avail))) avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(matrix(unlist(avail), ncol = length(avail))) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    P[!rows] <- 1
    return(P)
  }
  
  # ############################################################ #
  #### functionality="estimate/prediction/conditionals/raw" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","prediction","conditionals","raw")){
    
    # Store useful values
    if(functionality=="raw" && length(choiceVar)==1 && is.na(choiceVar)) choiceVar = alternatives[1]
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    nestnames <- names(nlStructure)
    if(!("root" %in% names(nlNests))) nlNests["root"] <- 1
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }
    
    ### Reorder V and avail to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(!anyNA(avail)) if(any(altnames != names(avail))) avail <- avail[altnames]
    
    combined_elements="root"
    j=1
    while(j<= length(nlStructure)){
      combined_elements=c(combined_elements,nlStructure[[j]])
      j=j+1
    }
    
    # Order tree structure
    nlStructure_ordered=list()
    element_list="root"
    j=1
    while(j>0){
      k=1
      temp=rep(TRUE,length(element_list))
      while(k<= length(element_list)){
        if(element_list[k] %in% altnames) temp[k]=FALSE
        k=k+1
      }
      element_list=element_list[temp]
      if(length(element_list)>0){
        nlStructure_ordered[[element_list[1]]]=nlStructure[[element_list[1]]]
        element_list=c(element_list,nlStructure[[element_list[1]]])
        element_list=element_list[-1]
      }
      j=length(element_list)
    }
    nlStructure=nlStructure_ordered
    
    # extract chosen utility
    chosenV <- Reduce('+',
                      lapply(as.list(1:nAlts),
                             FUN=function(i) (choiceVar==altcodes[i])*V[[altnames[i]]])
    )
    
    # substract chosen utility from all others for numerical stability
    V = lapply(X=V, FUN=function(v) v-chosenV)
    
    # loop over nests to create new utility elements and new availability terms
    k=length(nlStructure)
    while(k>0){
      nestK <- names(nlStructure)[k]
      V[[nestK]] = 0
      avail[[nestK]] = 1*( Reduce('+', avail[ nlStructure[[k]] ]) > 0 ) # calculate availability of nest
      j = 1
      while(j<= length(nlStructure[[k]])){
        nodeJ <- nlStructure[[k]][j]
        V[[nestK]] = V[[nestK]] + avail[[nodeJ]]*exp( V[[nodeJ]]/nlNests[[nestK]] )
        j = j+1
      }
      V[[nestK]] = nlNests[[nestK]]*log(V[[nestK]])
      k = k-1
    }
    
    # work out ancestors
    ancestors=list()
    j=1
    while(j<= length(altnames)){
      altJ <- altnames[[j]]
      ancestors[[altJ]] = altJ
      current = altJ
      k = length(nlStructure)
      while(k>0){
        if(current %in% nlStructure[[k]]){
          ancestors[[altnames[[j]]]] = c(ancestors[[altJ]],names(nlStructure)[k])
          current = names(nlStructure)[k]
        }
        k=k-1
      }
      j=j+1
    }
    
    # calculate log(probabilities)
    logPalts=list()
    j=1
    while(j <= length(altnames)){
      logPalts[[j]]=0
      k=1
      ancestorsJ <- ancestors[[altnames[[j]]]]
      while(k< length(ancestorsJ)){ # loop to level just below root
        current_V = V[[ ancestorsJ[k] ]]
        next_V    = V[[ ancestorsJ[k+1] ]]
        logPalts[[j]] = logPalts[[j]] + (current_V-next_V)/nlNests[[ ancestorsJ[k+1] ]]
        k=k+1
      }
      j=j+1
    }
    
    Palts = lapply(X=logPalts, FUN=exp)
    
    tempnames=names(V)
    names(Palts)=tempnames[1:length(altnames)]
    
    # consider availabilities (it assumes Palts and avail are in the same order)
    Palts <- mapply('*', Palts, avail[1:length(altnames)], SIMPLIFY = FALSE)
    
    if(functionality=="prediction"|(functionality=="raw")){
      Palts <- lapply(Palts, function(x) {x[is.na(x)] <- 0
      return(x)}) # replace all NaN by 0
      # add an additional column with chosen
      if(functionality=="prediction") Palts[["chosen"]] = Reduce('+', tmp <- mapply(function(x,y) (choiceVar==y)*x, Palts, as.list(altcodes), SIMPLIFY=FALSE) ) # keep only probability for chosen alternative
      Palts <- lapply(Palts, function(p) {
        # change likelihood of excluded columns for 1
        if(is.vector(p)) p[!rows]  <- ifelse(functionality=="prediction",NA,1)
        if(is.matrix(p)) p[!rows,] <- ifelse(functionality=="prediction",NA,1)
        if(is.array(p) & length(dim(p))==3) p[!rows,,] <- ifelse(functionality=="prediction",NA,1)
        return(p)
      })
    } else {
      Palts <- lapply(Palts, function(x) {x[is.na(x)] <- 0
      return(x)}) # replace all NaN by 0
      Palts = Reduce('+', tmp <- mapply(function(x,y) (choiceVar==y)*x, Palts, as.list(altcodes), SIMPLIFY=FALSE) ) # keep only probability for chosen alternative
      if(is.vector(Palts)) Palts[!rows]  <- 1
      if(is.matrix(Palts)) Palts[!rows,] <- 1
      if(is.array(Palts) && length(dim(Palts))==3) Palts[!rows,,] <- 1
    }
    
    return(Palts)
  }
  
  # ############################## #
  #### functionality="output" ####
  # ############################## #
  
  if(functionality=="output"){
    
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(V)
    avail_set <- FALSE
    root_set  <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    nestnames <- names(nlStructure)
    if(!("root" %in% names(nlNests))){
      root_set <- TRUE
      nlNests["root"] <- 1
    }
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }
    
    ### Reorder V and avail to match altnames order, if necessary
    if(any(altnames != names(V))) V <- V[altnames]
    if(!anyNA(avail)) if(any(altnames != names(avail))) avail <- avail[altnames]
    
    combined_elements="root"
    j=1
    while(j<= length(nlStructure)){
      combined_elements=c(combined_elements,nlStructure[[j]])
      j=j+1
    }
    
    # Order tree structure
    nlStructure_ordered=list()
    element_list="root"
    j=1
    while(j>0){
      k=1
      temp=rep(TRUE,length(element_list))
      while(k<= length(element_list)){
        if(element_list[k] %in% altnames) temp[k]=FALSE
        k=k+1
      }
      element_list=element_list[temp]
      if(length(element_list)>0){
        nlStructure_ordered[[element_list[1]]]=nlStructure[[element_list[1]]]
        element_list=c(element_list,nlStructure[[element_list[1]]])
        element_list=element_list[-1]
      }
      j=length(element_list)
    }
    nlStructure=nlStructure_ordered
    
    # extract chosen utility
    chosenV <- Reduce('+',
                      lapply(as.list(1:nAlts),
                             FUN=function(i) (choiceVar==altcodes[i])*V[[altnames[i]]])
    )
    
    # substract chosen utility from all others for numerical stability
    V = lapply(X=V, FUN=function(v) v-chosenV)
    
    # loop over nests to create new utility elements and new availability terms
    k=length(nlStructure)
    while(k>0){
      nestK <- names(nlStructure)[k]
      V[[nestK]] = 0
      avail[[nestK]] = 1*( Reduce('+', avail[ nlStructure[[k]] ]) > 0 ) # calculate availability of nest
      j = 1
      while(j<= length(nlStructure[[k]])){
        nodeJ <- nlStructure[[k]][j]
        V[[nestK]] = V[[nestK]] + avail[[nodeJ]]*exp( V[[nodeJ]]/nlNests[[nestK]] )
        j = j+1
      }
      V[[nestK]] = nlNests[[nestK]]*log(V[[nestK]])
      k = k-1
    }
    
    # work out ancestors
    ancestors=list()
    j=1
    while(j<= length(altnames)){
      altJ <- altnames[[j]]
      ancestors[[altJ]] = altJ
      current = altJ
      k = length(nlStructure)
      while(k>0){
        if(current %in% nlStructure[[k]]){
          ancestors[[altnames[[j]]]] = c(ancestors[[altJ]],names(nlStructure)[k])
          current = names(nlStructure)[k]
        }
        k=k-1
      }
      j=j+1
    }
    
    # calculate log(probabilities)
    logPalts=list()
    j=1
    while(j <= length(altnames)){
      logPalts[[j]]=0
      k=1
      ancestorsJ <- ancestors[[altnames[[j]]]]
      while(k< length(ancestorsJ)){ # loop to level just below root
        current_V = V[[ ancestorsJ[k] ]]
        next_V    = V[[ ancestorsJ[k+1] ]]
        logPalts[[j]] = logPalts[[j]] + (current_V-next_V)/nlNests[[ ancestorsJ[k+1] ]]
        k=k+1
      }
      j=j+1
    }
    
    Palts = lapply(X=logPalts, FUN=exp)
    
    tempnames=names(V)
    names(Palts)=tempnames[1:length(altnames)]
    
    # consider availabilities (it assumes Palts and avail are in the same order)
    Palts <- mapply('*', Palts, avail[1:length(altnames)], SIMPLIFY = FALSE)
    
    Palts <- lapply(Palts, function(x) {x[is.na(x)] <- 0
    return(x)}) # replace all NaN by 0
    P = Reduce('+', tmp <- mapply(function(x,y) (choiceVar==y)*x, Palts, as.list(altcodes), SIMPLIFY=FALSE) ) # keep only probability for chosen alternative
    if(is.vector(P)) P[!rows]  <- 1
    if(is.matrix(P)) P[!rows,] <- 1
    if(is.array(P) && length(dim(P))==3) P[!rows,,,drop=FALSE] <- 1
    
    # Construct summary table of availabilities and market share
    avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # n times each alt is available
    choicematrix = matrix(0,nrow=4,ncol=length(altnames))
    choicematrix[1,] = availprint
    j=1
    while(j<= length(altnames)){
      choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows) # times chosen
      j=j+1
    }
    choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
    choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
    choicematrix[4,!is.finite(choicematrix[4,])] <- 0
    rownames(choicematrix) = c("Times available","Times chosen","Percentage chosen overall","Percentage chosen when available")
    colnames(choicematrix) = altnames
    
    # Load function to write tree structure
    print_tree=function(nlStructure, ancestors){
      
      print_tree_level = function(nlStructure, component, preceding_nest_layer, space){
        j=1
        if(preceding_nest_layer!=0) space=c(space,"  |")
        while(j<=length(nlStructure[[component]])){
          space <- gsub("[']", " ", space)
          if(j==length(nlStructure[[component]])) space[length(space)] <- gsub("[|]", "'", space[length(space)])
          if(nlStructure[[component]][j] %in% altnames){
            depth <- length(space)
            cat("\n",space,rep("-",3*(maxDepth-depth)),"-Alternative: ",nlStructure[[component]][j], sep="")
          } else {
            cat("\n",space,"-Nest: ",nlStructure[[component]][j]," (",round(nlNests[[nlStructure[[component]][j]]],4), ")", sep="")
            print_tree_level(nlStructure, nlStructure[[component]][j], preceding_nest_layer+1, space)
          }
          j=j+1
        }
      }
      
      maxDepth <- max(sapply(ancestors, length))-1
      cat("Nest: ",names(nlStructure)[[1]]," (",round(nlNests[[names(nlStructure)[[1]]]],4),")", sep="")
      
      print_tree_level(nlStructure, "root", preceding_nest_layer=0, space="|")
    }
    
    ## write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
    #apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
    #                            error=function(e){
    #                              cat("apollo_nl could not retrieve apollo_control. No diagnostics in output.\n")
    #                              return(NA)
    #                            } )
    #if(!(length(apollo_control)==1 && is.na(apollo_control))){
    #  fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
    #  fileName <- file.path(tempdir(),fileName)
    #  fileConn <- tryCatch( file(fileName, open="at"),
    #                        error=function(e){
    #                          cat('apollo_nl could not write diagnostics to temporary file. No diagnostics in output.\n')
    #                          return(NA)
    #                        })
    #  if(!anyNA(fileConn)){
    #    sink(fileConn)
    #    on.exit({if(sink.number()>0) sink(); close(fileConn)})
    #    if(apollo_control$noDiagnostics==FALSE){
    #      cat('Overview of choices for NL model component:\n')
    #      print(round(choicematrix,0))
    #      cat("\n")
    #    }
    #    cat("Structure of nested logit component:")
    #    print_tree(nlStructure, ancestors)
    #    cat('\n')
    #    if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
    #  }
    #}
    
    #content <- list(round(choicematrix,2),
    #                ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
    #                ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""),
    #                "Final nesting structure:",
    #                capture.output(print_tree(nlStructure, ancestors)))
    content <- list(round(choicematrix,2))
    if(any(choicematrix[4,]==0)) content[[length(content)+1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(choicematrix[4,]==1)) content[[length(content)+1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    if(root_set==TRUE) content[[length(content)+1]] <- "Notice: Root lambda parameter set to 1."
    content[[length(content)+1]] <- "Nested structure:"
    content[[length(content)+1]] <- capture.output(print_tree(nlStructure, ancestors))
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for NL model component:", content, apolloLog)
    
  }
  
  return(P)
}