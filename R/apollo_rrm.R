#' Calculates Random Regret Minimisation model probabilities
#'
#' Calculates the probabilities of a Random Regret Minimisation model and can also perform other operations based on the value of the \code{functionality} argument.
#'
#' @param rrm_settings List of inputs of the RRM model. It should contain the following.
#'                     \itemize{
#'                       \item \strong{\code{alternatives}}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                       \item \strong{\code{avail}}: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1. These can be scalars or vectors (of length equal to rows in the database). A user can also specify \code{avail=1} to indicate universal availability, or omit the setting completely.
#'                       \item \strong{\code{choiceVar}}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                       \item \strong{\code{rum_inputs}}: Named list of (optional) deterministic utilities. Utilities of the alternatives to be included in combined RUM-RRM models. Names of elements must match those in \code{alternatives.}
#'                       \item \strong{\code{regret_inputs}}: Named list of regret functions. This should contain one list per attribute, where these lists themselves contain two vectors, namely a vector of attributes (at the alternative level) and parameters (either generic or attribute specific). Zeros can be used for omitted attributes for some alternatives. The order for each attribute needs to be the same as the order in \code{alternatives.}.
#'                       \item \strong{\code{regret_scale}}: Named list of regret scales. This should have the same length as 'rrm_settings$regret_inputs' or be a single entry in the case of a generic scale parameter across regret attributes.
#'                       \item \strong{\code{choiceset_scaling}}: Vector. One entry per row in the database, often set to 2 divided by the number of available alternatives.
#'                       \item \strong{\code{rows}}: Boolean vector. Consideration of which rows to include. Length equal to the number of observations (nObs), with entries equal to TRUE for rows to include, and FALSE for rows to exclude. Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
#'                     }
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
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"components"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"gradient"}}: List containing the likelihood and gradient of the model component.
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the probability of the chosen alternative.
#'           \item \strong{\code{"preprocess"}}: Returns a list with pre-processed inputs, based on \code{rrm_settings}.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'           \item \strong{\code{"report"}}: Choice overview
#'           \item \strong{\code{"shares_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when only constants are estimated.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'         }
#' @export
#' @importFrom utils capture.output
apollo_rrm <- function(rrm_settings, functionality){
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE, silent=FALSE, analyticGrad=TRUE)) ))
  
  ### Set or extract componentName
  modelType   = "RRM"
  if(is.null(rrm_settings[["componentName"]])){
    rrm_settings[["componentName"]] = ifelse(!is.null(rrm_settings[['componentName2']]),
                                             rrm_settings[['componentName2']], modelType)
    test <- functionality=="validate" && rrm_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType, ' without a componentName.', 
                                 ' The name was set to "', rrm_settings[["componentName"]], '" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, rrm_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("SPECIFICATION ISSUE - Duplicated componentName found (", rrm_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  if( !is.null(apollo_inputs[[paste0(rrm_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load rrm_settings from apollo_inputs
    rrm_settings <- apollo_inputs[[paste0(rrm_settings$componentName, "_settings")]]
  } else { 
    ### Do pre-processing
    if(is.null(names(rrm_settings))) stop('SYNTAX ISSUE - All elements inside the inputs lists for model components must be named, e.g. rrm_settings=list(alternatives=c(...), avail=...).')
    if(is.null(rrm_settings[["componentName"]])) stop('SYNTAX ISSUE - The settings of at least one model component is missing the mandatory "componentName" object.')
    # functionality
    test <- functionality %in% c("estimate","prediction","validate","zero_LL","shares_LL","conditionals","output","raw","preprocess", "components", "gradient", "report")
    if(!test) stop("SYNTAX ISSUE - Non-permissable setting for \"functionality\" for model component \"",rrm_settings$componentName,"\"")
    
    # Check for mandatory inputs
    mandatory <- c("alternatives", "choiceVar", "regret_inputs")
    for(i in mandatory) if(!(i %in% names(rrm_settings))) stop('SYNTAX ISSUE - The inputs list for model component "', rrm_settings$componentName, 
                                                               '" needs to include an object called "', i,'"!')
    # Check for optional inputs (avail and rows)
    if(is.null(rrm_settings[["rows"]])) rrm_settings[["rows"]]="all"
    if(is.null(rrm_settings[['avail']])){
      rrm_settings[['avail']]=1
      if(!apollo_inputs$silent && functionality=='validate') apollo_print('Setting "avail" is missing, so full availability is assumed.')
    }
    
    ### Store useful values
    rrm_settings$altnames = names(rrm_settings$alternatives)
    rrm_settings$altcodes = rrm_settings$alternatives
    rrm_settings$nAlt     = length(rrm_settings$alternatives)
    rrm_settings$nObs <- tryCatch(if(!is.null(apollo_inputs$database)) nrow(apollo_inputs$database) else stop('x'),
                                  error=function(e){
                                    lenV <- sapply(rrm_settings$regrets, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                                    lenA <- sapply(rrm_settings$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
                                    lenC <- ifelse(!is.null(rrm_settings$choiceVar),length(rrm_settings$choiceVar),length(rrm_settings$choiceShares[[1]]))
                                    return(max(lenV, lenA, lenC))
                                  })
    
    ### Format checks
    # alternatives
    test <- is.vector(rrm_settings$alternatives) & !is.null(names(rrm_settings$alternatives))
    if(!test) stop("SYNTAX ISSUE - The \"alternatives\" argument for model component \"",rrm_settings$componentName,"\" needs to be a named vector")
    # avail
    test <- is.list(rrm_settings$avail) || (length(rrm_settings$avail)==1 && rrm_settings$avail==1)
    if(!test) stop("SYNTAX ISSUE - The \"avail\" argument for model component \"",rrm_settings$componentName,"\" needs to be a list or set to 1")
    if(is.list(rrm_settings$avail)){
      lenA <- sapply(rrm_settings$avail, function(v) ifelse(is.array(v), dim(v)[1], length(v)) )
      test <- all(lenA==rrm_settings$nObs | lenA==1)
      if(!test) stop("SYNTAX ISSUE - All entries in \"avail\" for model component \"",rrm_settings$componentName,"\" need to be a scalar or a vector with one entry per observation in the \"database\"")
    }
    # choiceVar
    test <- is.vector(rrm_settings$choiceVar) && (length(rrm_settings$choiceVar)==rrm_settings$nObs || length(rrm_settings$choiceVar)==1)
    if(!test) stop("SYNTAX ISSUE - The \"choiceVar\" argument for model component \"",rrm_settings$componentName,"\" needs to be a scalar or a vector with one entry per observation in the \"database\"")
    
    # rows
    test <- is.vector(rrm_settings$rows)
    test <- test && ( is.logical(rrm_settings$rows) || all(rrm_settings$rows %in% 0:1) )
    test <- test && length(rrm_settings$rows)==rrm_settings$nObs
    test <- test || ( is.character(rrm_settings$rows) && length(rrm_settings$rows)==1 && rrm_settings$rows=="all" )
    if(!test) stop("SYNTAX ISSUE - The 'rows' argument for model component '", rrm_settings$componentName, 
                   "' needs to be \"all\" or a vector of logical statements or 0/1 entries", 
                   " with one entry per observation in the 'database'")
    if( all(rrm_settings$rows %in% c(0,1)) ) (rrm_settings$rows <- rrm_settings$rows>0)
    
    ### Expand availabilities if necessary
    rrm_settings$avail_set <- FALSE
    if(length(rrm_settings$avail)==1 && rrm_settings$avail==1){
      rrm_settings$avail <- as.list(setNames(rep(1,rrm_settings$nAlt), rrm_settings$altnames))
      rrm_settings$avail_set <- TRUE
    }
    
    ### Check that avail is available for all alternatives
    if(!all(rrm_settings$altnames %in% names(rrm_settings$avail))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"", rrm_settings$componentName,"\" do not match those in \"avail\".")
    
    ### Reorder availabilities
    rrm_settings$avail <- rrm_settings$avail[rrm_settings$altnames]
    
    
    ### Check RRM specific inputs
    rrm_settings$altnames      = names(rrm_settings$alternatives)
    ### Check if all rum & regret inputs are character (text strings) or not
    isText <- is.null(rrm_settings$rum_inputs) || is.list(rrm_settings$rum_inputs)
    isText <- isText && all(sapply(rrm_settings$rum_inputs, is.character))
    isText <- isText && is.list(rrm_settings$regret_inputs)
    for(le in rrm_settings$regret_inputs){
      isText <- isText && is.list(le) && length(le)==2 && !is.null(names(le))
      isText <- isText && all(names(le) %in% c("x", "b"))
      isText <- isText && all(sapply(le[["x"]], is.character))
      isText <- isText && all(sapply(le[["b"]], is.character))
    }
    isText <- isText && (is.null(rrm_settings$regret_inputs) || 
                           (is.list(rrm_settings$regret_scale) && all(sapply(rrm_settings$regret_scale, is.character))))

    ### If rum & regret inputs are not all text, make them text and run pre-processing
    if(!isText){
      ap <- tryCatch(get("apollo_probabilities", envir=parent.frame(), inherits=TRUE),
                     error=function(e) NULL)
      ab <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                     error=function(e) NULL)
      af <- tryCatch(get("apollo_fixed", envir=parent.frame(), inherits=TRUE),
                     error=function(e) NULL)
      if(is.null(ap) || is.null(ab) || is.null(af)) stop("INTERNAL ISSUE - Could not fetch apollo_probabilities", 
                                                         " or apollo_beta when pre-processing", 
                                                         " RRM component")
      tmpSilent <- apollo_inputs$silent
      tmpAnGrad <- apollo_inputs$apollo_control$analyticGrad
      apollo_inputs$silent <- TRUE
      apollo_inputs$apollo_control$analyticGrad <- FALSE
      L <- apollo_modifyUserDefFunc(ab, af, ap, apollo_inputs, validate=FALSE)
      apollo_inputs$apollo_lcPars    <- L$apollo_lcPars
      apollo_inputs$apollo_randCoeff <- L$apollo_randCoeff
      apollo_inputs$apollo_scaling   <- L$apollo_scaling
      ap <- L$apollo_probabilities
      newSet <- ap(ab, apollo_inputs, functionality="preprocess")
      test <- is.list(newSet) && !is.null(names(newSet))
      test <- test && (paste0(rrm_settings[["componentName"]], "_settings") %in% names(newSet))
      if(!test) stop("INTERNAL ISSUE - Could not pre-process RRM component")
      newSet <- newSet[[paste0(rrm_settings[["componentName"]], "_settings")]]
      if(functionality=="preprocess") return(newSet)
      apollo_inputs$silent <- tmpSilent
      apollo_inputs$apollo_control$analyticGrad <- tmpAnGrad
      rrm_settings <- newSet
      rm(ap, ab, af, tmpSilent, tmpAnGrad, L, newSet, test)
    } else {
      ### create empty RUM part if not provided
      if(is.null(rrm_settings$rum_inputs)){
        rrm_settings$rum_inputs        = as.list(rep("0",length(rrm_settings$altnames)))
        names(rrm_settings$rum_inputs) = rrm_settings$altnames  
      }
      ### check names
      if(any(!(names(rrm_settings$rum_inputs)%in%rrm_settings$altnames))) stop("INTERNAL ISSUE - Some inputs in 'rrm_settings$rum_inputs' use names for alternatives not defined in rrm_settings$alternatives!")
      ### add missing parts to RUM part
      if(length(rrm_settings$rum_inputs)<length(rrm_settings$altnames)){
        missing_alts = subset(rrm_settings$altnames,!(rrm_settings$altnames%in%names(rrm_settings$rum_inputs)))
        for(j in missing_alts) rrm_settings$rum_inputs[[j]]="0"}
      ### check format of regret_inputs
      if(!is.list(rrm_settings$regret_inputs)) stop("SYNTAX ISSUE - The element 'rrm_settings$regret_inputs' needs to be a list!")
      for(s in 1:length(rrm_settings$regret_inputs)){
        ### check that all entries in regret_inputs are lists with x and b
        if(!(length(rrm_settings$regret_inputs[[s]])&&(all(c("x","b")%in%names(rrm_settings$regret_inputs[[s]]))))) stop("SYNTAX ISSUE - The element ",names(rrm_settings$regret_inputs)[s]," in 'rrm_settings$regret_inputs' needs to be a list with two elements, named 'x' and 'b'!")
        ### check that length of all x entries is equal to the number of alternatives
        if(!(length(rrm_settings$regret_inputs[[s]][["x"]])==length(rrm_settings$altnames))) stop("SYNTAX ISSUE - The element ",names(rrm_settings$regret_inputs)[s],"$x in 'rrm_settings$regret_inputs' needs to have the same length as the number of alternatives in your model. If attribute ",names(rrm_settings$regret_inputs)[s]," is not used for a specific alternative, you should use \"0\" for the corresponding entry!")
        ### check that length of all b entries is 1 or equal to the number of alternatives
        if(!(length(rrm_settings$regret_inputs[[s]][["b"]])%in%c(1,length(rrm_settings$altnames)))) stop("SYNTAX ISSUE - The element ",names(rrm_settings$regret_inputs)[s],"$b in 'rrm_settings$regret_inputs' needs to either have the same length as the number of alternatives in your model or be a single entry in the case of generic parameters!")
        ### expand all b entries to the be the same length as x
        if(length(rrm_settings$regret_inputs[[s]][["b"]])==1) rrm_settings$regret_inputs[[s]][["b"]]=rep(rrm_settings$regret_inputs[[s]][["b"]],length(rrm_settings$altnames))
      }
      ### check for missing scale part
      if(is.null(rrm_settings$regret_scale)) rrm_settings$regret_scale="1"
      ### check for length of scale part
      if(!(length(rrm_settings$regret_scale)%in%c(1,length(rrm_settings$regret_inputs)))) stop("SYNTAX ISSUE - The element regret_scale in 'rrm_settings' needs to either have the same length as 'rrm_settings$regret_inputs' or be a single entry in the case of a generic scale parameter across regret attributes!")
      ### expand scale to be the same length as regret_inputs
      if(length(rrm_settings$regret_scale)==1) rrm_settings$regret_scale=rep(rrm_settings$regret_scale,length(rrm_settings$regret_inputs))
      
      ### 25 Jan added choiceset_scaling
      ### check for missing choiceset_scaling
      if(is.null(rrm_settings$choiceset_scaling)) rrm_settings$choiceset_scaling="1"
      
      ### initialise regret function by starting with the RUM part (which may be just be zeros)
      rrm_settings$regrets = rrm_settings$rum_inputs
      ### loop over alternatives
      for(i in 1:length(rrm_settings$altnames)){
        ### loop over attributes to add regret part
        for(k in 1:length(rrm_settings$regret_inputs)){
          ### loop over other alternatives  
          J=1:length(rrm_settings$altnames)
          J=subset(J,!(J%in%i))
          for(j in J){
            ### availabilities need to be used twice, once inside the log to avoid numerical issues depending on coding of unavailable options, and once outside to ensure that unavailable options don't contribute to regret
            ### 25 Jan 2023 - added avail for i inside the comparison
            ### 25 Jan 2023 - added choiceset scaling
            rrm_settings$regrets[[i]]=paste0(rrm_settings$regrets[[i]]," - ",rrm_settings$choiceset_scaling," * ",rrm_settings$regret_scale[[k]]," * rrm_settings$avail[[\"",rrm_settings$altnames[j],"\"]]* log(1+exp(1/",rrm_settings$regret_scale[[k]],"*(rrm_settings$avail[[\"",rrm_settings$altnames[j],"\"]]*(",rrm_settings$regret_inputs[[k]]$b[[j]],"*",rrm_settings$regret_inputs[[k]]$x[[j]],")-rrm_settings$avail[[\"",rrm_settings$altnames[i],"\"]]*(",rrm_settings$regret_inputs[[k]]$b[[i]],"*",rrm_settings$regret_inputs[[k]]$x[[i]],"))))")
            #rrm_settings$regrets[[i]]=paste0(rrm_settings$regrets[[i]]," + ",rrm_settings$regret_scale[k]," * log(1+exp(1/",rrm_settings$regret_scale[k],"*(rrm_settings$avail[[\"",rrm_settings$altnames[j],"\"]]*",rrm_settings$regret_inputs[[k]]$b[j],"*",rrm_settings$regret_inputs[[k]]$x[j],"-",rrm_settings$regret_inputs[[k]]$b[i],"*",rrm_settings$regret_inputs[[k]]$x[i],")))")
          }
        }
        ### add function() at the beginning
        rrm_settings$regrets[[i]]=paste0("function() (",rrm_settings$regrets[[i]],")")
      }
      
      # Turn text into functions and insert scaling in them
      rrm_settings$regrets = lapply(rrm_settings$regrets,str2lang)
      rrm_settings$regrets = lapply(rrm_settings$regrets,eval)
      
      
      ### Expand rows if necessary, and update nObs
      if(length(rrm_settings$rows)==1 && rrm_settings$rows=="all") rrm_settings$rows <- rep(TRUE, rrm_settings$nObs)
      rrm_settings$nObs <- sum(rrm_settings$rows)
      # Filter rows in avail and choiceVar
      if(any(!rrm_settings$rows)){
        rrm_settings$avail <- lapply(rrm_settings$avail, 
                                     function(av) if(length(av)==1) return(av) else return(av[rrm_settings$rows]))
        rrm_settings$choiceVar <- apollo_keepRows(rrm_settings$choiceVar, rrm_settings$rows)
      }
      
      ### Create Y
      rrm_settings$Y <- lapply(as.list(rrm_settings$alternatives), function(i) rrm_settings$choiceVar==i)
      
      # Record availability of chosen alternative
      rrm_settings$chosenAvail <- Reduce('+', mapply('*', rrm_settings$Y, rrm_settings$avail, SIMPLIFY=FALSE))
      
      # Determine which rrm likelihood to use (R or C++)
      if(apollo_inputs$apollo_control$cpp & !apollo_inputs$silent) apollo_print("No C++ optimisation for RRM available")
      # Using R likelihood
      rrm_settings$probs_RRM=function(rrm_settings, all=FALSE){
        # Fix choiceVar if "raw" and choiceVar==NA
        rrm_settings$choiceNA = FALSE
        if(all(is.na(rrm_settings$choiceVar))){
          rrm_settings$choiceVar = rrm_settings$alternatives[1]
          rrm_settings$choiceNA = TRUE
        }
        # Set regrets of unavailable alternatives to 0 to avoid numerical issues (eg attributes = -999)
        rrm_settings$regrets <- mapply(function(r,a) apollo_setRows(r, !a, 0), 
                                       rrm_settings$regrets, rrm_settings$avail, SIMPLIFY=FALSE)
        # if probabilities for all alternatives are requested, then P is a list
        if(all){
          ### work with subtracting the chosenR
          chosenR <- mapply("*", rrm_settings$Y, rrm_settings$regrets, SIMPLIFY=FALSE)
          chosenR <- Reduce('+', chosenR)
          rrm_settings$regrets <- lapply(X=rrm_settings$regrets, "-", chosenR)
          rrm_settings$regrets <- lapply(X=rrm_settings$regrets, FUN=exp)
          # consider availabilities (it assumes V and avail are in the same order)
          rrm_settings$regrets <- mapply('*', rrm_settings$regrets, rrm_settings$avail, SIMPLIFY = FALSE)
          # calculate the denominator of the Logit probability expression for RRM
          denom = Reduce('+',rrm_settings$regrets)
          P <- lapply(rrm_settings$regrets, "/", denom)
          if(any(sapply(P, anyNA))){
            P <- lapply(P, function(p){p[is.na(p)] <- 1; return(p)} )
            P <- mapply('*', P, rrm_settings$avail, SIMPLIFY = FALSE)
            sP <- Reduce("+", P)
            P <- lapply(P, "/", sP)
          }
          if(!rrm_settings$choiceNA) P[["chosen"]] <- Reduce("+", mapply("*", rrm_settings$Y, P, SIMPLIFY=FALSE))
          # if only the probability of the chosen alternative is requested, then P is vector or a 3-dim array
        } else { 
          ### work with subtracting the chosenR
          chosenR <- mapply("*", rrm_settings$Y, rrm_settings$regrets, SIMPLIFY=FALSE)
          chosenR <- Reduce('+', chosenR)
          rrm_settings$regrets <- lapply(X=rrm_settings$regrets, "-", chosenR)
          rrm_settings$regrets <- lapply(X=rrm_settings$regrets, FUN=exp)
          # consider availabilities (it assumes regrets and avail are in the same order)
          rrm_settings$regrets <- mapply('*', rrm_settings$regrets, rrm_settings$avail, SIMPLIFY = FALSE)
          # calculate the denominator of the Logit probability expression
          denom = Reduce('+',rrm_settings$regrets)
          P <- rrm_settings$chosenAvail/denom
        }
        return(P)
      }
      
      # Create diagnostics function for RRM
      rrm_settings$rrm_diagnostics <- function(inputs, apollo_inputs, data=TRUE, param=TRUE){
        
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
          #if(inputs$avail_set) apollo_print("Availability not provided (or some elements are NA). Full availability assumed.", type="i")
          apollo_print("\n")
          apollo_print(paste0('Overview of choices for RRM model component ', 
                              ifelse(inputs$componentName=='model', '', inputs$componentName), ':'))
          print(round(choicematrix,2))
        }
        
        return(invisible(TRUE))
      }
      
      # Store model type
      rrm_settings$modelType <- modelType
      
      # Construct necessary input for gradient (including gradient of utilities)
      apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                              error=function(e) return(NULL))
      test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
      test <- test && all(sapply(rrm_settings$regrets, is.function))
      test <- test && apollo_inputs$apollo_control$analyticGrad
      rrm_settings$gradient <- FALSE
      if(test){
        rrm_settings$regrets  <- rrm_settings$regrets[rrm_settings$altnames] # reorder regrets
        rrm_settings$dR       <- apollo_dVdB(apollo_beta, apollo_inputs, rrm_settings$regrets)
        rrm_settings$gradient <- !is.null(rrm_settings$dR)
      }; rm(test)
      
      # Return rrm_settings if pre-processing
      if(functionality=="preprocess"){
        # Remove things that change from one iteration to the next
        # Nothing changes
        return(rrm_settings)
      }
    }
  }
  
  # ################################################## #
  #### Transform regrets into numeric and drop rows ####
  # ################################################## #

  ### Execute regrets (makes sure we are now working with vectors/matrices/arrays and not functions)
  if(any(sapply(rrm_settings$regrets, is.function))){
    avail_backup=rrm_settings$avail
    rrm_settings$avail=lapply(rrm_settings$avail,apollo_insertRows,r=rrm_settings$rows,val=1)
    e <- environment()
    rrm_settings$regrets = lapply(rrm_settings$regrets, function(f){
      if(is.function(f)){ environment(f) <- e; return(f()) } else return(f)
    } )
    rm(e)
    rrm_settings$avail=avail_backup
    rm(avail_backup)
    #cat("\n dim(regrets)=", paste0(sapply(rrm_settings$regrets, length), collapse=", "), "\n", sep="")
  } 
  rrm_settings$regrets <- lapply(rrm_settings$regrets, function(r) if(is.matrix(r) && ncol(r)==1) as.vector(r) else r)
  
  ### Reorder regrets and drop rows if necessary
  rrm_settings$regrets <- rrm_settings$regrets[rrm_settings$altnames]
  if(!all(rrm_settings$rows)) rrm_settings$regrets <- lapply(rrm_settings$regrets, apollo_keepRows, r=rrm_settings$rows)

  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if(functionality=="validate"){
    # Check that alternatives are named in altcodes and regrets
    if(is.null(rrm_settings$altnames) || is.null(rrm_settings$altcodes) || is.null(names(rrm_settings$regrets))) stop("SYNTAX ISSUE - Alternatives for model component \"",rrm_settings$componentName,"\" must be named, both in 'alternatives' and 'regrets'.")
    
    if(!apollo_inputs$apollo_control$noValidation){
      # Check there are no repeated alternatives names
      if(length(unique(rrm_settings$altnames))!=length(rrm_settings$altnames)) stop('SYNTAX ISSUE - Names of alternatives must be unique. Check definition of "alternatives".')
      
      # Check that there are at least two alternatives
      minAlts <- 2
      if(rrm_settings$nAlt<minAlts) stop("SYNTAX ISSUE - Model component \"",rrm_settings$componentName,"\"  requires at least ", minAlts, " alternatives")
      
      # Check that choice vector is not empty
      if(length(rrm_settings$choiceVar)==0) stop("SYNTAX ISSUE - Choice vector is empty for model component \"",rrm_settings$componentName,"\"")
      
      if(rrm_settings$nObs==0) stop("SYNTAX ISSUE - No data for model component \"",rrm_settings$componentName,"\"")
      
      # Check regrets and avail elements are named correctly
      if(!all(rrm_settings$altnames %in% names(rrm_settings$regrets))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"",rrm_settings$componentName,"\" do not match those in \"regrets\".")
      if(!all(rrm_settings$altnames %in% names(rrm_settings$avail))) stop("SYNTAX ISSUE - The names of the alternatives for model component \"",rrm_settings$componentName,"\" do not match those in \"avail\".")
      
      # Check that there are no values in the choice column for undefined alternatives
      rrm_settings$choiceLabs <- unique(rrm_settings$choiceVar)
      if(!all(rrm_settings$choiceLabs %in% rrm_settings$altcodes)) stop("SYNTAX ISSUE - The data contains values for \"choiceVar\" for model component \"",rrm_settings$componentName,"\" that are not included in \"alternatives\".")
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(rrm_settings$avail)) if( !all(unique(rrm_settings$avail[[i]]) %in% 0:1) ) stop("SYNTAX ISSUE - Some availability values for model component \"",rrm_settings$componentName,"\" are not 0 or 1.")
      # check that at least 2 alternatives are available in at least one observation
      if(max(Reduce('+',rrm_settings$avail))==1) stop("CALCULATION ISSUE - Only one alternative is available for each observation for model component \"",rrm_settings$componentName,"!")
      # check that nothing unavailable is chosen
      for(j in 1:rrm_settings$nAlt) if(any(rrm_settings$choiceVar==rrm_settings$altcodes[j] & rrm_settings$avail[[j]]==0)){
        txt <-  paste0('WARNING: The data contains cases where alternative ', rrm_settings$altnames[j], 
                       ' is chosen for model component "',rrm_settings$componentName, '" despite being', 
                       ' listed as unavailable. This will cause the chosen probability to be', 
                       ' zero, and potentially lead to an invalid LL.')
        apollo_print(txt)
      }
      
      # Check that no available alternative has regret = NA
      # Requires setting non available alternatives utility to 0 first
      rrm_settings$regrets <- mapply(function(r,a) apollo_setRows(r, !a, 0), rrm_settings$regrets, rrm_settings$avail, SIMPLIFY=FALSE)
      if(!all(sapply(rrm_settings$regrets, function(r) all(is.finite(r))))) stop('CALCULATION ISSUE - Some regrets for model component "',
                                                                     rrm_settings$componentName, 
                                                                     '" contain values that are not finite numbers!')
    } 
    
    if(!apollo_inputs$apollo_control$noDiagnostics) rrm_settings$rrm_diagnostics(rrm_settings, apollo_inputs)
    
    testL = rrm_settings$probs_RRM(rrm_settings, all=FALSE)
    if(any(!rrm_settings$rows)) testL <- apollo_insertRows(testL, rrm_settings$rows, 1)
    if(all(testL==0)) stop("CALCULATION ISSUE - All observations have zero probability at starting value for model component \"",rrm_settings$componentName,"\"")
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0("Some observations have zero probability at starting value for model component \"",rrm_settings$componentName,"\"", sep=""), type="i")
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:length(rrm_settings$avail)) if(length(rrm_settings$avail[[i]])==1) rrm_settings$avail[[i]] <- rep(rrm_settings$avail[[i]], rrm_settings$nObs)
    nAvAlt <- rowSums(do.call(cbind, rrm_settings$avail)) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    if(any(!rrm_settings$rows)) P <- apollo_insertRows(P, rrm_settings$rows, 1)
    return(P)
  }

  # ############################### #
  #### functionality="shares_LL" ####
  # ############################### #
  
  if(functionality=="shares_LL"){
    for(i in 1:length(rrm_settings$avail)) if(length(rrm_settings$avail[[i]])==1) rrm_settings$avail[[i]] <- rep(rrm_settings$avail[[i]], rrm_settings$nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(do.call(cbind, rrm_settings$avail)) # number of available alts in each observation
    Y = do.call(cbind,rrm_settings$Y)
    if(var(nAvAlt)==0){
      Yshares = colSums(Y)/nrow(Y)
      P = as.vector(Y%*%Yshares)
    } else {
      ## Estimate model with constants only
      rrm_ll = function(b, A, Y) as.vector(Y%*%c(b,0) - log(rowSums( A%*%exp(c(b,0)) )))
      A = do.call(cbind, rrm_settings$avail)
      b = maxLik::maxLik(rrm_ll, start=rep(0, rrm_settings$nAlt - 1), 
                         method='BFGS', finalHessian=FALSE, A=A, Y=Y)$estimate
      P = exp(rrm_ll(b, A, Y))
    }
    if(any(!rrm_settings$rows)) P <- apollo_insertRows(P, rrm_settings$rows, 1)
    return(P)
  }
  
  # ################################################################# #
  #### functionality="estimate/prediction/conditionals/raw/output" ####
  # ################################################################# #
  
  if(functionality %in% c("estimate","conditionals", "components", "output")){
    P <- rrm_settings$probs_RRM(rrm_settings, all=FALSE)
    if(any(!rrm_settings$rows)) P <- apollo_insertRows(P, rrm_settings$rows, 1)
    return(P)
  }
  
  if(functionality %in% c("prediction","raw")){
    P <- rrm_settings$probs_RRM(rrm_settings, all=TRUE)
    if(any(!rrm_settings$rows)) P <- lapply(P, apollo_insertRows, r=rrm_settings$r, val=NA)
    return(P)
  }
  
  # ############################## #
  #### functionality="gradient" ####
  # ############################## #
  
  if(functionality=="gradient"){
    # Verify everything necessary is available
    if(is.null(rrm_settings$dR) || !all(sapply(rrm_settings$dR, is.function))) stop("INTERNAL ISSUE - Analytical gradient could not be calculated. Please set apollo_control$analyticGrad=FALSE.")
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) stop("INTERNAL ISSUE - apollo_rrm could not fetch apollo_beta for gradient estimation."))
    if(is.null(apollo_inputs$database)) stop("INTERNAL ISSUE - apollo_rrm could not fetch apollo_inputs$database for gradient estimation.")
    
    # Calculate probabilities and derivatives of utilities for all alternatives
    P    <- rrm_settings$probs_RRM(rrm_settings, all=TRUE)
    Pcho <- P[["chosen"]]
    P    <- P[-which(names(P)=="chosen")]
    e <- list2env(c(as.list(apollo_beta), apollo_inputs$database, list(apollo_inputs=apollo_inputs), rrm_settings), hash=TRUE)
    for(i in 1:length(rrm_settings$dR)) environment(rrm_settings$dR[[i]]) <- e
    dR<- lapply(rrm_settings$dR, function(dr) dr()) # One element per alt. Each element is an nPar long list
    if(!all(rrm_settings$rows)) for(i in 1:length(dR)) dR[[i]] <- lapply(dR[[i]], apollo_keepRows, rrm_settings$rows)
    for(i in 1:rrm_settings$nAlt) dR[[i]] <- lapply(dR[[i]], 
                                                    function(drik){ # Make dR=0 for unavailable alternatives
                                                      test <- length(drik)==1 && length(rrm_settings$avail[[i]])>1
                                                      if(test) drik <- rep(drik, rrm_settings$nObs)
                                                      drik[!rrm_settings$avail[[i]]] <- 0
                                                      return(drik)
                                                    })
    
    # Calculate gradient
    GA<- mapply(function(y,p) y - p, rrm_settings$Y, P, SIMPLIFY=FALSE)
    G <- list()
    for(k in 1:length(dR[[1]])){
      dRk   <- lapply(dR, function(dr) dr[[k]])
      G[[k]] <- Reduce("+", mapply("*", GA, dRk, SIMPLIFY=FALSE))
    }
    G <- lapply(G, "*", Pcho)
    
    # Restore rows and return
    if(!all(rrm_settings$rows)){
      Pcho <- apollo_insertRows(Pcho, rrm_settings$rows, 1)
      G    <- lapply(G, apollo_insertRows, r=rrm_settings$rows, val=0)
    }
    return(list(like=Pcho, grad=G))
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(rrm_settings$rrm_diagnostics(rrm_settings, apollo_inputs, param=FALSE))
    P$param <- capture.output(rrm_settings$rrm_diagnostics(rrm_settings, apollo_inputs, data =FALSE))
    return(P)
  }
  
}
