#' Calculates MDCEV likelihoods.
#'
#' Calculates the likelihood of a Multiple Discrete Continuous Extreme Value (MDCEV) model.
#'
#' @param mdcev_settings List of settings for the MDCEV model. It must include the following.
#'                       \itemize{
#'                         \item \strong{\code{V}}: Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#'                         \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'V'.
#'                         \item \strong{\code{alpha}}: Named list. Alpha parameters for each alternative, including for any outside good. As many elements as alternatives.
#'                         \item \strong{\code{gamma}}: Named list. Gamma parameters for each alternative, excluding any outside good. As many elements as inside good alternatives.
#'                         \item \strong{\code{sigma}}: Numeric scalar. Scale parameter of the model extreme value type I error.
#'                         \item \strong{\code{cost}}: Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{avail}}: Named list. Availabilities of alternatives, one element per alternative. Names of elements must match those in argument 'alternatives'. Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation.
#'                         \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{budget}}: Numeric vector. Budget for each observation.
#'                         \item \strong{\code{minConsumption}}: Named list of scalars or numeric vectors. Minimum consumption of the alternatives, if consumed. As many elements as alternatives. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{outside}}: Character. Optional name of the outside good.
#'                         \item \strong{\code{rows}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                         \item \strong{\code{componentName}}: Character. Name given to model component.
#'                       }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"} Used for model estimation.
#'                        \item \code{"prediction"} Used for model predictions.
#'                        \item \code{"validate"} Used for validating input.
#'                        \item \code{"zero_LL"} Used for calculating null likelihood.
#'                        \item \code{"conditionals"} Used for calculating conditionals.
#'                        \item \code{"output"} Used for preparing output after model estimation.
#'                        \item \code{"raw"} Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the observed consumption for each observation.
#'           \item \strong{\code{"prediction"}}: A matrix with one row per observation, and columns indicating means and s.d. of continuous and discrete predicted consumptions.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'         }
#' @export
apollo_mdcev <- function(mdcev_settings,functionality){
  if(is.null(mdcev_settings[["componentName"]])       ) mdcev_settings[["componentName"]]="MDCEV"
  componentName    = mdcev_settings[["componentName"]]
  
    # Check values in mdcev_settings
  if( is.null(mdcev_settings[["alternatives"]])     ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"alternatives\"!")
  if( is.null(mdcev_settings[["avail"]])            ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"avail\"!")
  if( is.null(mdcev_settings[["continuousChoice"]]) ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"continuousChoice\"!")
  if( is.null(mdcev_settings[["V"]])                ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"V\"!")
  if( is.null(mdcev_settings[["alpha"]])            ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"alpha\"!")
  if( is.null(mdcev_settings[["gamma"]])            ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"gamma\"!")
  if( is.null(mdcev_settings[["sigma"]])            ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"sigma\"!")
  if( is.null(mdcev_settings[["cost"]])             ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"cost\"!")
  if( is.null(mdcev_settings[["budget"]])           ) stop("The mdcev_settings list for model component \"",componentName,"\" needs to include an object called \"budget\"!")
  if( is.null(mdcev_settings[["minConsumption"]])   ) mdcev_settings[["minConsumption"]] = NA
  if( is.null(mdcev_settings[["rows"]])             ) mdcev_settings[["rows"]]           = "all"
  if( is.null(mdcev_settings[["outside"]])          ) mdcev_settings[["outside"]]        = NA
  
  # Extract values in mdcev_settings
  alternatives     = mdcev_settings[["alternatives"]]
  avail            = mdcev_settings[["avail"]]
  continuousChoice = mdcev_settings[["continuousChoice"]]
  V                = mdcev_settings[["V"]]
  alpha            = mdcev_settings[["alpha"]]
  gamma            = mdcev_settings[["gamma"]]
  sigma            = mdcev_settings[["sigma"]]
  cost             = mdcev_settings[["cost"]]
  budget           = mdcev_settings[["budget"]]
  minConsumption   = mdcev_settings[["minConsumption"]]
  rows             = mdcev_settings[["rows"]]
  outsideName      = mdcev_settings[["outside"]]
  
  if(anyNA(outsideName) && "outside" %in% alternatives) outsideName = "outside"
  if(!anyNA(outsideName) && !(outsideName %in% alternatives)) stop("Name provided for outside good for model component \"",componentName,"\" does not correspond to an alternative!")
  
  # # Basic checks and useful variables
  if(!is.list(continuousChoice) || !all(sapply(continuousChoice, is.numeric)) ) stop("'continuousChoice' for model component \"",componentName,"\" should be a list of numerical vectors, each with as many rows as observations.")
  nObs  <- max(sapply(continuousChoice, length))
  nAlts <- length(V)
  if(!is.vector(alternatives) || !is.character(alternatives)) stop("'alternatives' for model component \"",componentName,"\" must be a character vector.")
  if(length(alternatives)!=nAlts) stop("The length of 'alternatives' for model component \"",componentName,"\" does not match with the length of 'V'.")
  if(!is.list(avail) || any(sapply(avail, function(x) !is.logical(x) & !is.numeric(x) & !(length(x) %in% c(1,nObs)) )) ) stop("'avail' for model component \"",componentName,"\" should be a list of scalars or vectors (each with as many elements as observations), with values 0 or 1.")
  if(!is.list(V) || !all(sapply(V, is.numeric)) ) stop("'V' for model component \"",componentName,"\" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.")
  if(!is.list(alpha) || !all(sapply(alpha, is.numeric)) ) stop("'alpha' for model component \"",componentName,"\" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.")
  if(!is.list(gamma) || !all(sapply(gamma, is.numeric)) ) stop("'gamma' for model component \"",componentName,"\" should be a list of numerical vectors/matrices/arrays, each with as many rows as observations.")
  if(!is.list(cost) || !all(sapply(cost, is.numeric)) ) stop("'cost' for model component \"",componentName,"\" should be a list of numerical scalar or vectors (each with as many rows as observations), with positive values.")
  if(!is.numeric(budget) || !(length(budget) %in% c(1,nObs)) ) stop("'budget' for model component \"",componentName,"\" should be a scalar or a vector with as many elements as observations.")
  if(!(length(rows)==1 && rows=="all") && !(is.vector(rows) && is.logical(rows) && length(rows)==nAlts)) stop("The argument \"rows\" for model component \"",componentName,"\" needs to either be \"all\" or a boolean vector of length equal to the number of rows in the data!")
  if(!anyNA(minConsumption)) if(!is.list(minConsumption) || length(minConsumption)!=nObs || !all(sapply(minConsumption, is.numeric)) || !all(sapply(minConsumption, length) %in% c(1,nObs))){
    stop("'minConsumption' for model component \"",componentName,"\", if included, must be a list of numeric scalars or vectors (each with as many elelemnts as observations).")
  }
  
  
  # Necesary for apollo_mdcevOutside and apollo_mdcevInside to "see" apollo_inputs
  apollo_inputs <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE ), error=function(e) return(NA))
  
  # Add outside good's gamma if omitted
  if(!anyNA(outsideName) && is.null(gamma[[outsideName]])) gamma[[outsideName]]=1
  
  # Put the outsidegood first, if it exists
  if(!anyNA(outsideName)) alternatives=c(outsideName,alternatives[-which(alternatives==outsideName)])
  
  # Sort lists based on the order of 'alternatives'
  if(any(alternatives != names(V))) V <- V[alternatives]
  if(any(alternatives != names(avail))) avail <- avail[alternatives]
  if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
  if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
  if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
  if(any(alternatives != names(cost))) cost <- cost[alternatives]
  if(!anyNA(minConsumption)) if(any(alternatives != names(minConsumption))) minConsumption <- minConsumption[alternatives]
  
  # If there is an outside good, rename it to "outside"
  if(!anyNA(outsideName)){
    alternatives = c("outside", alternatives[2:length(alternatives)])
    names(V    ) = alternatives
    names(avail) = alternatives
    names(alpha) = alternatives
    names(gamma) = alternatives
    names(cost ) = alternatives
    names(continuousChoice) = alternatives
    if(!anyNA(minConsumption)) names(minConsumption) = alternatives 
  }
  
  # Expand "rows"
  if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs)
  
  # For excluded rows: put all consumption in first alternative, and set budget and cost to 1 (unless scalar)
  # No need no change budget and cost if they are scalar, as invalid number will fail for considered rows too.
  if(length(budget   )==nObs){ budget[!rows]    <- 1; budgetTmp <- 1 } else budgetTmp <- budget   
  if(length(cost[[1]])==nObs){ cost[[1]][!rows] <- 1; costTmp   <- 1 } else costTmp   <- cost[[1]]
  continuousChoice[[1]][!rows] <- budgetTmp/costTmp
  for(j in 2:nAlts){
    continuousChoice[[j]][!rows] <- 0
    if(length(cost[[j]])==nObs) cost[[j]][!rows] <- 1
  } 
  
  # For unavailable alternatives: Set consumption to zero and cost to 1 (unless scalar).
  # No need no change cost if scalar, as an invalid number will fail for wvailable rows too.
  for(j in 1:nAlts) if(length(avail[[j]])==nObs){
    continuousChoice[[j]][!avail[[j]]] <- 0
    if(length(cost[[j]])==nObs) cost[[j]][!avail[[j]]] <- 1
  }
  
  # set minConsumption to NA if they are all zero
  if(!anyNA(minConsumption)) if(Reduce('+',lapply(minConsumption, sum))==0) minConsumption=NA
  
  # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
  V <- mapply(function(v,a) apollo_setRows(v, !a | !rows, 0), V, avail, SIMPLIFY=FALSE)
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  if(functionality=="validate"){
    # Create availability if needed
    avail_set <- FALSE
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(TRUE,nAlts)), alternatives)
    }
    
    # Useful variables
    discrete_choice <- lapply(continuousChoice, function(x) x>0)
    if(!anyNA(apollo_inputs)) apollo_control <- apollo_inputs$apollo_control else{
      apollo_control <- list(noValidation=FALSE, noDiagnostics=FALSE)
    }
    
    ### VALIDATION
    if(apollo_control$noValidation==FALSE){
      
      # Check that sigma is not random (actually, it could be, but it leads to weird results)
      if(!is.vector(sigma)) stop("Sigma for model component \"",componentName,"\" should not be random")
      if(length(sigma)!=1 && length(sigma)!=nObs) stop("Sigma for model component \"",componentName,"\" should be either a scalar or a vector with as many elements as observations")
      
      # Check there are at least two alternatives
      if(nAlts<2) stop("Model component \"",componentName,"\" requires at least two products")
      
      # Check that choice vector is not empty
      if(nObs==0) stop("No data for model component \"",componentName,"\"")
      
      # Check that first product is outside good
      if(!anyNA(outsideName) && alternatives[1]!="outside") stop("First product for model component \"",componentName,"\" must be called \"outside\"!")
      
      namesinside=names(V)[names(V)!="outside"]
      # Check labels
      if(!all(alternatives %in% names(V))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"V\"!")
      if(!all(alternatives %in% names(alpha))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"alpha\"!")
      if(!all(namesinside %in% names(gamma))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"gamma\"!")
      if(!all(alternatives %in% names(continuousChoice))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"continuousChoice\"!")
      if(!all(alternatives %in% names(cost))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"cost\"!")
      if(!all(alternatives %in% names(avail))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"avail\"!")
      
      # check that nothing unavailable is chosen
      for(j in 1:length(V)) if( any( discrete_choice[[j]] & avail[[j]]==0 ) ) stop("Product", j, "chosen despite being listed as unavailable for model component \"",componentName,"\"!")
      
      # check that outside good is always chosen
      if(!anyNA(outsideName)){
        txt <- paste0("Outside good ", ifelse(outsideName!="outside", paste0("(",outsideName,") "), ""),"for model component \"",componentName,"\" should always be chosen!")
        if(any(continuousChoice[["outside"]]<=0)) stop(txt)
      }
      
      # check that all costs are positive
      if( sum(sapply(cost, function(x) sum(x<=0))) > 0 ) stop("Costs for model component \"",componentName,"\" must be strictly positive for all products!")
      
      # check consumption is non negative for all products
      for(i in 1:length(continuousChoice)) if(any(continuousChoice[[i]]<0)){
        stop( paste0("Consumption values of alternative", alternatives[i], "for model component \"",componentName,"\" must be non-negative!") )
      }
      
      # Check budget>0
      if( (length(budget)==1 && budget<=0) || (length(budget)>1 && any(budget[rows]<=0)) ) stop("Budget for model component \"",componentName,"\" for some rows in data is less than or equal to zero!")
      
      # check that full budget is consumed in each row, nothing more, nothing less
      expenditure <- Reduce("+", mapply("*", continuousChoice, cost, SIMPLIFY=FALSE))
      if(any(abs(expenditure-budget)>10^-10)) stop("Expenditure for some observations for model component \"",componentName,"\" is either less or more than budget!")
      
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1 for model component \"",componentName,"\".")
      
      # check that availability of outside is always 1
      if(!anyNA(outsideName) && any(avail[[1]]!=1) ) stop("Outside good is not available for some observations for model component \"",componentName,"\". It should always be available.")
      
      # check that if minimum consumption exists, it has the same names as alternatives, and that no consumptions are less than minConsumption if alternative is available
      if(!anyNA(minConsumption)){
        if(!all(alternatives %in% names(minConsumption))) stop("Labels in \"alternatives\" for model component \"",componentName,"\" do not match those in \"minConsumption\"!")
        for(i in 1:length(alternatives)) if(any(continuousChoice[[i]][avail[[i]] & rows] < minConsumption[[i]][avail[[i]] & rows])){
          stop( paste0("Consumption of alternative ", alternatives[i], " for model component \"",componentName,"\" is smaller than its listed minConsumption") )
        }
      }
      
      # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
      V <- mapply(function(v,a) apollo_setRows(v, !a | !rows, 0), V, avail, SIMPLIFY=FALSE)
      if(any(sapply(V, anyNA))) cat("\nThe utility of at least one alternative for model component \"",componentName,"\" contains one or more NA values")
    }
    
    if(apollo_control$noDiagnostics==FALSE){  
      availprint = colSums(matrix(unlist(avail), ncol = length(avail))[rows,])
      
      choicematrix = matrix(0,nrow=4,ncol=length(V))
      
      choicematrix[1,] = availprint
      for(j in 1:length(V)){
        choicematrix[2,j]=sum(discrete_choice[[j]][rows]==1)
        choicematrix[3,j]=sum(continuousChoice[[j]][rows])/choicematrix[1,j]
        choicematrix[4,j]=sum(continuousChoice[[j]][rows])/choicematrix[2,j]
      }
      choicematrix <- apply(choicematrix, MARGIN=2, function(x) {x[!is.finite(x)] <- 0; return(x)})
      rownames(choicematrix) = c("Times available","Observations in which chosen",
                                 "Average consumption when available",
                                 "Average consumption when chosen")
      if(!anyNA(outsideName)) colnames(choicematrix) = c(outsideName, names(V)[2:length(V)]) else{
        colnames(choicematrix) = alternatives
      }
      
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog(title   = paste0("Overview of choices  for model component \"",componentName,"\":"), 
                    content = content, apolloLog)
      
    }
    
    testL=apollo_mdcev(mdcev_settings, functionality="estimate")
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"",componentName,"\"")
    if(any(testL==0)) cat("\nSome observations have zero probability at starting value for model component \"",componentName,"\"")
    return(invisible(testL))
  }
  
  
  # ################################# #
  #### functionality != "validate" ####
  # ################################# #
  # Call apollo_mdcevOutside or apollo_mdcevInside
  if(!anyNA(outsideName)){
    ans <- apollo_mdcevOutside(V, alternatives, alpha, gamma, sigma, cost, 
                               avail, continuousChoice, budget, functionality, 
                               minConsumption, outsideName, rows, componentName)
  } else {
    ans <- apollo_mdcevInside(V, alternatives, alpha, gamma, sigma, cost, 
                              avail, continuousChoice, budget, functionality, 
                              minConsumption, rows, componentName)
  }
  
  return(ans)
}
