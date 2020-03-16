#'
#' Calculates probabilities of a multinomial logit model.
#'
#' @param mnl_settings List of inputs of the MNL model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                       \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{\code{V}}: Named list of deterministic utilities . Utilities of the alternatives. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component.
#'                     }
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
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @examples
#' ### Load data
#' data(apollo_modeChoiceData)
#' database <- apollo_modeChoiceData
#' rm(apollo_modeChoiceData)
#'
#' ### Parameters
#' b = list(asc_1=0, asc_2=0, asc_3=0, asc_4=0, tt=0, tc=0, acc=0)
#'
#' ### List of utilities
#' V = list()
#' V[['car' ]] = b$asc_1 + b$tt*database$time_car  + b$tc*database$cost_car
#' V[['bus' ]] = b$asc_2 + b$tt*database$time_bus  + b$tc*database$cost_bus  +
#'               b$acc*database$access_bus
#' V[['air' ]] = b$asc_3 + b$tt*database$time_air  + b$tc*database$cost_air  +
#'               b$acc*database$access_air
#' V[['rail']] = b$asc_4 + b$tt*database$time_rail + b$tc*database$cost_rail +
#'               b$acc*database$access_rail
#'
#' ### MNL settings
#' mnl_settings <- list(
#'    alternatives = c(car=1, bus=2, air=3, rail=4),
#'    avail        = list(car=database$av_car, bus=database$av_bus,
#'                        air=database$av_air, rail=database$av_rail),
#'    choiceVar    = database$choice,
#'    V            = V
#' )
#'
#' ### Compute choice probabilities using MNL model
#' apollo_mnl(mnl_settings, functionality="estimate")
#' @importFrom Deriv Deriv
apollo_mnl_2 <- function(mnl_settings, functionality){
  if(is.null(mnl_settings[["alternatives"]])) stop("The mnl_settings list needs to include an object called \"alternatives\"!")
  if(is.null(mnl_settings[["avail"]])) stop("The mnl_settings list needs to include an object called \"avail\"!")
  if(is.null(mnl_settings[["choiceVar"]])) stop("The mnl_settings list needs to include an object called \"choiceVar\"!")
  if(is.null(mnl_settings[["V"]])) stop("The mnl_settings list needs to include an object called \"V\"!")
  if(is.null(mnl_settings[["rows"]])) mnl_settings[["rows"]]="all"
  
  ### Useful values
  alternatives = mnl_settings[["alternatives"]]
  avail        = mnl_settings[["avail"]]
  choiceVar    = mnl_settings[["choiceVar"]]
  V            = mnl_settings[["V"]]
  rows         = mnl_settings[["rows"]]
  altnames     = names(alternatives)
  altcodes     = alternatives
  nAlt         = length(alternatives)

  
  ### if V is a function, copy function into a new object, and create V as values
  
  ### in what follows, next steps only happen if previous is true
  # 1: test if any elements of V are a function, if not, do nothing
  test=sapply(V,is.function)
  if(any(test)){
    # 2: if any elements are a function, test that all are a function, else fail
    if(any(!test)) stop("Some elements in V are functions, while others are not!")
    # 3: if functionality==gradient, create a copy of V in Vfunction
    if(functionality=="gradient") Vfunction=V
    # 4: replace each element in the V list by the output of the relevant function
    V=lapply(V,function(v) v())
  }

  nObs <- tryCatch(nrow( get("apollo_inputs", parent.frame(), inherits=FALSE)$database ),
                   error=function(e){
                     lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenA <- sapply(avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                     lenC <- length(choiceVar)
                     return(max(lenV, lenA, lenC))
                   })
  
  ### Format checks
  # alternatives
  test <- is.vector(alternatives) & !is.null(names(alternatives))
  if(!test) stop("The \"alternatives\" argument needs to be a named vector")
  # avail
  test <- is.list(avail) || (length(avail)==1 && avail==1)
  if(!test) stop("The \"avail\" argument needs to be a list or set to 1")
  if(is.list(avail)){
    lenA <- sapply(avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
    test <- all(lenA==nObs | lenA==1)
    if(!test) stop("All entries in \"avail\" need to be a scalar or a vector with one entry per observation in the \"database\"")
  }
  # choiceVar
  test <- is.vector(choiceVar) && (length(choiceVar)==nObs || length(choiceVar)==1)
  if(!test) stop("The \"choiceVar\" argument needs to be a scalar or a vector with one entry per observation in the \"database\"")
  # V
  if(!is.list(V)) stop("The \"V\" argument needs to be a list")
  lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
  test <- all(lenV==nObs | lenV==1)
  if(!test) stop("Each element of \"V\" needs to be a scalar or a vector/matrix/cube with one row per observation in the \"database\"")  
  # rows
  test <- is.vector(rows) && ( (is.logical(rows) && length(rows)==nObs) || (length(rows)==1 && rows=="all") )
  if(!test) stop("The \"rows\" argument needs to be \"all\" or a vector of boolean statements with one entry per observation in the \"database\"")
  # functionality
  test <- functionality %in% c("estimate","prediction","gradient","validate","zero_LL","conditionals","output","raw")
  if(!test) stop("Non-permissable setting for \"functionality\"")
  
  ### Expand availabilities if =1
  avail_set <- FALSE
  if(length(avail)==1 && avail==1){
    avail <- as.list(setNames(rep(1,nAlt), altnames))
    avail_set <- TRUE
  } 
  
  ### Filter rows
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, length(nObs))
  if(any(!rows)){
    avail <- lapply(avail, function(av) if(length(av)==1) return(av) else return(av[rows]))
    choiceVar <- choiceVar[rows]
    V <- lapply(V, apollo_keepRows, r=rows)
    # update nObs
    lenV <- sapply(V, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
    lenA <- sapply(avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
    lenC <- length(choiceVar)
    nObs <- max(lenV, lenA, lenC)
    rm(lenV, lenA, lenC)
  }
  
  ### Reorder arguments
  avail <- avail[altnames]
  V     <- V[altnames]
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Store useful values
    apollo_control <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_control,
                               error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    
    # check rows statement
    #if(length(rows)!=length(choiceVar)) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of the rows in the data!")
    
    # Check that alternatives are named in altcodes and V
    if(is.null(altnames) || is.null(altcodes) || is.null(names(V))) stop("Alternatives must be named, both in 'alternatives' and 'V'.")
    
    if(apollo_control$noValidation==FALSE){
      # Check that there are at least two alternatives
      if(nAlt<2) stop("MNL requires at least two alternatives")
      
      # Check that choice vector is not empty
      if(length(choiceVar)==0) stop("No choices to model")
      if(nObs==0) stop("No utilities defined")
      
      choiceLabs <- unique(choiceVar)
      if(!all(altnames %in% names(V))) stop("Alternatives' names in \"alternatives\" do not match those in \"V\".")
      if(!all(altnames %in% names(avail))) stop("Alternatives' names in \"alternatives\" do not match those in \"avail\".")
      
      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop('A value in "choiceVar" column is not included in "alternatives".')
      
      # check that nothing unavailable is chosen
      for(j in 1:nAlt) if(any(choiceVar==altcodes[j] & avail[[j]]==0)) stop("Alternative ",altnames[j]," chosen despite being listed as unavailable\n")
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")
      
      # Check that no available alternative has utility = NA
      # Requires setting non available alternatives utility to 0 first
      V <- mapply(function(v,a) apollo_setRows(v, !a, 0), V, avail, SIMPLIFY=FALSE)
      if(any(sapply(V, anyNA))) cat("\nAt least one utility contains one or more NA values at parameter starting values")
    }
    
    if(apollo_control$noDiagnostics==FALSE){
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # Construct summary table of availabilities and market share
      choicematrix = matrix(0, nrow=4, ncol=length(altnames), 
                            dimnames=list(c("Times available", "Times chosen", "Percentage chosen overall",
                                            "Percentage chosen when available"), altnames))
      choicematrix[1,] = unlist(lapply(avail, sum))
      for(j in 1:length(altnames)) choicematrix[2,j] = sum(choiceVar==altcodes[j]) # number of times each alt is chosen
      choicematrix[3,] = choicematrix[2,]/nObs*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      choicematrix[4,!is.finite(choicematrix[4,])] <- 0
      
      content <- list(round(choicematrix,2),
                      ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
                      ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""))
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided to 'apollo_mnl' (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for MNL model component:", content, apolloLog)
    }
    
    return(invisible(TRUE))
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #

  if(functionality=="gradient"){
    apollo_beta = get("apollo_beta" , parent.frame(), inherits=TRUE)
    apollo_fixed= get("apollo_fixed", parent.frame(), inherits=TRUE)
    freeparams  = names(apollo_beta[!names(apollo_beta)%in%apollo_fixed])
    K <- length(freeparams)
    J <- length(Vfunction)
    #for(i in 1:J) environment(Vfunction[[i]]) <- baseenv()
    
    dV<- lapply(Vfunction, Deriv::Deriv, x=freeparams)
    dV<- lapply(dV, function(f){
      nLi <- length(body(f)) # the search could be made recursive
      if( all( freeparams %in% names(body(f)) ) & as.character(body(f)[[1]])=="c" ){
        body(f)[[1]] <- substitute(cbind)
      } else body(f)[[nLi]][[1]] <- substitute(cbind)
      return(f)
    })
    #for(i in 1:J) environment(dV[[i]]) <- baseenv()
    
    G = matrix(0, nrow=nObs, ncol=K, dimnames=list(NULL, freeparams))
    Y = lapply(as.list(alternatives), function(i) choiceVar==i)
    g <- function(b, db){
      # Change functions' environment to a new environment with b and db
      e <- list2env(c(as.list(b), db), hash=TRUE)
      for(i in 1:length(dV)){
        environment(Vfunction[[i]]) <- e
        environment(dV[[i]])        <- e
      }
      # Calculate probabilities
      V <- lapply(Vfunction, function(v) v())
      chosenV <- Reduce("+", mapply("*", Y, V, SIMPLIFY=FALSE))
      V <- lapply(X=V, FUN=function(v) exp(v-chosenV)) # avoids numerical issues
      V <- mapply("*", V, avail, SIMPLIFY=FALSE)
      P <- lapply(V, "/", Reduce("+", V))
      # Calculate gradient (won't work with mixing)
      GA<- mapply(function(y,p) y - p, Y, P, SIMPLIFY=FALSE)
      dV<- lapply(dV, function(dv) dv())
      for(k in 1:K){
        dVk   <- lapply(dV, function(dv) dv[,k])
        G[,k] <- Reduce("+", mapply(function(ga,dvk) ga*dvk, GA, dVk, SIMPLIFY=FALSE))
      }
      return(G)
    }
    environment(g) <- new.env(parent=baseenv())
    assign( "G",  G, envir=environment(g))
    assign( "Y",  Y, envir=environment(g))
    assign( "K",  K, envir=environment(g))
    assign("dV", dV, envir=environment(g))
    assign("avail", avail, envir=environment(g))
    assign("Vfunction", Vfunction, envir=environment(g))
    
    return(g)
  }
  
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, avail)) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    if(any(!rows)) P <- apollo_insertRows(P, rows, 1)
    return(P)
  }
  
  # ############################################################ #
  #### functionality="estimate/prediction/conditionals/raw" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate","prediction","conditionals","raw")){
    
    # Fix choiceVar if "raw" and choiceVar==NA
    choiceNA = FALSE
    if(all(is.na(choiceVar))){
      choiceVar = alternatives[1]
      choiceNA = TRUE
    }
    
    # Set utility of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
    V <- mapply(function(v,a) apollo_setRows(v, !a, 0), V, avail, SIMPLIFY=FALSE)
    
    # Extract chosen utility
    chosenV <- Reduce('+',
                      lapply(as.list(1:nAlt),
                             FUN=function(i) (choiceVar==altcodes[i])*V[[altnames[i]]])
    )
    
    # Extract chosen avail, needed for choices of unavailable
    chosenAvail <- Reduce('+',
                          lapply(as.list(1:nAlt),
                                 FUN=function(i) (choiceVar==altcodes[i])*avail[[altnames[i]]])
    )
    
    # substract chosen utility from all others for numerical stability
    V = lapply(X=V, FUN=function(v) v-chosenV)
    
    # exponentiate utilities
    V = lapply(X=V, FUN=exp)
    
    # consider availabilities (it assumes V and avail are in the same order)
    V <- mapply('*', V, avail, SIMPLIFY = FALSE)
    
    # calculate the denominator of the logit probability expression
    # now simply the addition of all exponentiated utilities
    denom = Reduce('+',V)
    
    # calculate probabilities
    # if probabilities for all alternatives are requested, then P is a list
    # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
    if((functionality=="prediction")|(functionality=="raw")){
      P <- lapply(V, function(ev) ev/denom)
      if(!choiceNA) P[["chosen"]]= chosenAvail/denom
      if(any(!rows)) P <- lapply(P, apollo_insertRows, r=rows, val=NA) # insert excluded rows with value 1
    } else {
      P <- chosenAvail/denom
      if(any(!rows)) P <- apollo_insertRows(P, rows, 1) # insert excluded rows with value 1
    }
    
    return(P)
  }
  
  # ############################## #
  #### functionality="output" ####
  # ############################## #
  
  if(functionality=="output"){
    P <- apollo_mnl_2(mnl_settings, functionality="estimate")
    
    # turn scalar availabilities into vectors
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
    
    # Construct summary table of availabilities and market share
    choicematrix = matrix(0, nrow=4, ncol=length(altnames), 
                          dimnames=list(c("Times available", "Times chosen", "Percentage chosen overall",
                                          "Percentage chosen when available"), altnames))
    choicematrix[1,] = unlist(lapply(avail, sum))
    for(j in 1:length(altnames)) choicematrix[2,j] = sum(choiceVar==altcodes[j]) # number of times each alt is chosen
    choicematrix[3,] = choicematrix[2,]/nObs*100 # market share
    choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
    choicematrix[4,!is.finite(choicematrix[4,])] <- 0
    
    content <- list(round(choicematrix,2),
                    ifelse(any(choicematrix[4,]==0), "Warning: some alternatives are never chosen in your data!", ""),
                    ifelse(any(choicematrix[4,]==1), "Warning: some alternatives are always chosen when available!", ""))
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided to 'apollo_mnl' (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for MNL model component:", content, apolloLog)
    
    return(P)
  }
}
