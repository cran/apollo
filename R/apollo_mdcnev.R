#' Calculates MDCNEV likelihoods with an outside good.
#' 
#' Calculates the likelihood of a Multiple Discrete Continuous Nested Extreme Value (MDCNEV) model with an outside good.
#' 
#' @param mdcnev_settings List of settings for the MDCEV model. It must include the following.
#'                       \itemize{
#'                         \item \strong{\code{V}}: Named list. Utilities of the alternatives. Names of elements must match those in argument 'alternatives'.
#'                         \item \strong{\code{alternatives}}: Character vector. Names of alternatives, elements must match the names in list 'V'.
#'                         \item \strong{\code{alpha}}: Named list. Alpha parameters for each alternative, including for the outside good. As many elements as alternatives.
#'                         \item \strong{\code{gamma}}: Named list. Gamma parameters for each alternative, including for the outside good. As many elements as alternatives.
#'                         \item \strong{\code{mdcnevNests}}: Named list. Lambda parameters for each nest. Elements must be named with the nest name. The lambda at the root is fixed to 1, and therefore must be no be defined. The value of the estimated mdcnevNests parameters should be between 0 and 1 to ensure consistency with random utility maximization.
#'                         \item \strong{\code{mdcnevStructure}}: Numeric matrix. One row per nest and one column per alternative. Each element of the matrix is 1 if an alternative belongs to the corresponding nest.
#'                         \item \strong{\code{cost}}: Named list of numeric vectors. Price of each alternative. One element per alternative, each one as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{avail}}: Named list. Availabilities of alternatives, one element per alternative. Names of elements must match those in argument 'alternatives'. Value for each element can be 1 (scalar if always available) or a vector with values 0 or 1 for each observation. If all alternatives are always available, then user can just omit this argument.
#'                         \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount of consumption of each alternative. One element per alternative, as long as the number of observations or a scalar. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{budget}}: Numeric vector. Budget for each observation.
#'                         \item \strong{\code{minConsumption}}: Named list of scalars or numeric vectors. Minimum consumption of the alternatives, if consumed. As many elements as alternatives. Names must match those in \code{alternatives}.
#'                         \item \strong{\code{outside}}: Character. Alternative name for the outside good. Default is "outside"
#'                         \item \strong{\code{rows}}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                         \item \strong{\code{componentName}}: Character. Name given to model component.
#'                       }
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
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the observed consumption for each observation.
#'           \item \strong{\code{"prediction"}}: A matrix with one row per observation, and columns indicating means and s.d. of continuous and discrete predicted consumptions.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}, but it also runs a set of tests to validate the function inputs.
#'           \item \strong{\code{"zero_LL"}}: Not implemented. Returns a vector of NA with as many elements as observations.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"estimate"}
#'         }
#' @export
#' @importFrom mnormt rmnorm
#' @importFrom stats setNames
apollo_mdcnev <- function(mdcnev_settings,functionality){
  if(is.null(mdcnev_settings[["componentName"]])) mdcnev_settings[["componentName"]]="MDCNEV"
  componentName     = mdcnev_settings[["componentName"]]
  
  if(is.null(mdcnev_settings[["alternatives"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"alternatives\"!")
  if(is.null(mdcnev_settings[["avail"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"avail\"!")
  if(is.null(mdcnev_settings[["continuousChoice"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"continuousChoice\"!")
  if(is.null(mdcnev_settings[["V"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"V\"!")
  if(is.null(mdcnev_settings[["alpha"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"alpha\"!")
  if(is.null(mdcnev_settings[["gamma"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"gamma\"!")
  if(is.null(mdcnev_settings[["mdcnevNests"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"mdcnevNests\"!")
  if(is.null(mdcnev_settings[["mdcnevStructure"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"mdcnevStructure\"!")
  if(is.null(mdcnev_settings[["cost"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"cost\"!")
  if(is.null(mdcnev_settings[["budget"]])) stop("The mdcnev_settings list for model component \"",componentName,"\" needs to include an object called \"budget\"!")
  if(is.null(mdcnev_settings[["minConsumption"]])) mdcnev_settings[["minConsumption"]]=NA
  if(is.null(mdcnev_settings[["rows"]])) mdcnev_settings[["rows"]]="all"
  if(is.null(mdcnev_settings[["outside"]])) mdcnev_settings[["outside"]]=NA ### new
  
  alternatives      = mdcnev_settings[["alternatives"]]
  avail             = mdcnev_settings[["avail"]]
  continuousChoice  = mdcnev_settings[["continuousChoice"]]
  V                 = mdcnev_settings[["V"]]
  alpha             = mdcnev_settings[["alpha"]]
  gamma             = mdcnev_settings[["gamma"]]
  sigma             = 1
  mdcnevNests       = mdcnev_settings[["mdcnevNests"]]
  mdcnevStructure   = mdcnev_settings[["mdcnevStructure"]]
  cost              = mdcnev_settings[["cost"]]
  budget            = mdcnev_settings[["budget"]]
  minConsumption    = mdcnev_settings[["minConsumption"]]
  rows              = mdcnev_settings[["rows"]]
  outsideName       = mdcnev_settings[["outside"]]   ### new
  
  if(anyNA(outsideName) && "outside" %in% alternatives) outsideName = "outside"  
  if(anyNA(outsideName)) stop("No outside good provided for MDCNEV model!")      
  if(is.null(gamma[[outsideName]])) gamma[[outsideName]]=1
  # Reorder 'alternatives'
  # original location of outside good
  outside_loc=which(alternatives==outsideName)
  alternatives=c(outsideName,alternatives[-outside_loc])
  
  # Reorder V and avail to match 'alternatives' order, if necessary
  if(any(alternatives != names(V))) V <- V[alternatives]
  if(any(alternatives != names(avail))) avail <- avail[alternatives]
  if(any(alternatives != names(alpha))) alpha <- alpha[alternatives]
  if(any(alternatives != names(gamma))) gamma <- gamma[alternatives]
  if(any(alternatives != names(continuousChoice))) continuousChoice <- continuousChoice[alternatives]
  if(any(alternatives != names(cost))) cost <- cost[alternatives]
  if(!anyNA(minConsumption)) if(any(alternatives != names(minConsumption))) minConsumption <- minConsumption[alternatives]
  mdcnevStructure=cbind(mdcnevStructure[,outside_loc],mdcnevStructure[,-outside_loc])
  
  # Rename first element
  alternatives[1] <- "outside"
  names(V    ) = alternatives
  names(avail) = alternatives
  names(alpha) = alternatives
  names(gamma) = alternatives
  names(cost ) = alternatives
  names(continuousChoice) = alternatives
  if(!anyNA(minConsumption)) names(minConsumption) = alternatives
  colnames(mdcnevStructure) = alternatives
  

  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){

    #if(!is.null(gamma[["outside"]])){
    #  cat("\nWarning: You provided a gamma for the outside good,\n         which will be ignored in the model.\n")
    #  gamma[["outside"]]=NULL
    #}

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

      # Check there are at least two alternatives
      if(nAlts<2) stop("MDCEV requires at least two products")

      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")

      # Check that first product is outside good
      if(alternatives[1]!="outside") stop("First product must be called \"outside\"!")

      namesinside=names(V)[names(V)!="outside"]
      
      # Check labels
      if(!all(alternatives %in% names(V))) stop("Labels in \"alternatives\" do not match those in \"V\"!")
      if(!all(alternatives %in% names(alpha))) stop("Labels in \"alternatives\" do not match those in \"alpha\"!")
      if(!all(namesinside %in% names(gamma))) stop("Labels in \"alternatives\" do not match those in \"gamma\"!")
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

      # check that outside good is always chosen
      if(sum(continuousChoice[[1]]==0)>0) stop("First product should always be chosen as it is an outside good!")

      # checks for consumption for individual products
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

      # check that full budget is consumed in each row, nothing more, nothing less
      if(sum(abs(budget_check-budget)>10^-10)) stop("Expenditure for some observations is either less or more than budget!")

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
      
      # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
      V <- mapply(function(v,a) apollo_setRows(v, !a | !rows, 0), V, avail, SIMPLIFY=FALSE)
      if(any(sapply(V, anyNA))) cat("\nAt least one utility contains one or more NA values")

      # checks that are specific to cnlStructure component
      if(nrow(mdcnevStructure)!=length(mdcnevNests)) stop("Tree structure needs one row per nest!")
      if(ncol(mdcnevStructure)!=nAlts) stop("Tree structure needs one column per alternative!")
      if(any(colSums(mdcnevStructure)!=1)) stop("Each alternative must be allocated to one nest only!")

      #cat("\nAll checks passed for MDCNEV model component\n")
    }

    if(apollo_control$noDiagnostics==FALSE){
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
      colnames(choicematrix) = c(outsideName, names(V)[2:length(V)])
      colnames(mdcnevStructure) <- c(outsideName, names(V)[2:length(V)])
      rownames(mdcnevStructure) <- names(mdcnevNests)
      
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                           "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      content[[length(content) + 1]] <- "Structure for MDCNEV model component:"
      content[[length(content) + 1]] <- mdcnevStructure
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for MDCEV model component:", content, apolloLog)

    }
    testL=apollo_mdcnev(mdcnev_settings, functionality="estimate")
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for model component \"",componentName,"\"")
    if(any(testL==0)) cat("\nSome observations have zero probability at starting value for model component \"",componentName,"\"")
    return(invisible(testL))
  }
  #if(!is.null(gamma[["outside"]])) gamma[["outside"]]=NULL
  
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
  # ############################################### #
  if(functionality %in% c("estimate", "conditionals", "raw")){
    omega <- mdcnevStructure
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    nNests<- length(mdcnevNests)
    avail_set <- FALSE
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      for(i in which(!rows)){
        continuousChoice[[j]][1] <- budget[i]
        for(j in 2:nAlts) continuousChoice[[j]][i] <- 0
      }
    }
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), alternatives)
    }

    
    # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
    V <- mapply(function(v,a) apollo_setRows(v, !a | !rows, 0), V, avail, SIMPLIFY=FALSE)
    #if(any(sapply(V, anyNA))) stop("At least one utility contains one or more NA values")

    term4base=function(t, ars, qrs){
      sumx=1
      if(qrs==0|qrs==1)
        sumx = 1
      if(qrs >1 & ars==1)
        sumx = 1
      if (ars==2)
        sumx = ((1-t)/t)*qrs*(qrs-1)/2
      if (qrs==3 & ars==3){
        sumx = (2*(1-t)/t)+1
        sumx = sumx * ((1-t)/t)}
      if (qrs==4 & ars==3){
        sumx = ((3*(1-t)/t)+1)*((2*(1-t)/t)+(1*(1-t)/t))
        sumx = sumx + (((2*(1-t)/t)+1)*((1-t)/t))}
      if (qrs==4 & ars==4)
        sumx = (((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t))
      if (qrs==5 & ars==3){
        sumx = ((4*(1-t)/t)+1) * ((3*(1-t)/t)+(2*(1-t)/t)+(1*(1-t)/t))
        sumx = sumx + ((3*(1-t)/t)+1) * ((2*(1-t)/t)+(1*(1-t)/t))
        sumx = sumx + (((2*(1-t)/t)+1)*(1*(1-t)/t))}
      if (qrs==5 & ars==4){
        sumx = ((4*(1-t)/t)+2)*((3*(1-t)/t)+1) * ((2*(1-t)/t)+(1*(1-t)/t))
        sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)
        sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)}
      if (qrs==5 & ars==5)
        sumx = ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*((1-t)/t)
      if (qrs==6 & ars==3){
        sumx = ((5*(1-t)/t)+1)*(10*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==6 & ars==4){
        sumx = ((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==6 & ars==5){
        sumx = ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==6 & ars==6)
        sumx = ((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
      if (qrs==7 & ars==3){
        sumx = ((6*(1-t)/t)+1)*(15*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+1)*(10*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==7 & ars==4){
        sumx = ((6*(1-t)/t)+2)*((5*(1-t)/t)+1)*(10*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(3*(1-t)/t)}
      if (qrs==7 & ars==5){
        sumx = ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((4*(1-t)/t)+1)*(6*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+3)*((5*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==7 & ars==6){
        sumx = ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((3*(1-t)/t)+1)*(3*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+4)*((5*(1-t)/t)+3)*((4*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((6*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
        sumx = sumx + ((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)}
      if (qrs==7 & ars==7) sumx = ((6*(1-t)/t)+5)*((5*(1-t)/t)+4)*((4*(1-t)/t)+3)*((3*(1-t)/t)+2)*((2*(1-t)/t)+1)*(1*(1-t)/t)
      return(sumx)
    }

    term4final=function(mdcnevNests, r_current_combo, q_person){
      out <- 1
      for(s in 1:length(r_current_combo)) if(q_person[s]>0) out <- out * term4base(mdcnevNests[[s]], r_current_combo[s], q_person[s])
      return(out)
    }

    term5=function(r_current_combo, q_person){
      out <- 0
      for(s in 1:length(r_current_combo)) if(q_person[s]>0) out <- out + q_person[s]- r_current_combo[s] + 1
      out <- factorial(out - 1)
      return(out)
    }

    
    
    discrete_choice=list()
    for(k in 1:nAlts) discrete_choice[[alternatives[k]]] = (continuousChoice[[k]]>0)*1

    q = list() # stores how many different products where purchased in each nest
    for(s in 1:nNests){
      alts   <- which(as.vector(omega[s,])>0)
      q[[s]] <- Reduce("+", discrete_choice[alts])
    }
    
    mdcnevNestsAlt = list() # called theta in old code
    for(k in 1:nAlts) mdcnevNestsAlt[[k]] <- mdcnevNests[[which(as.vector(omega[,k])>0)]]
    
    # Compute V
    V[[1]] = (alpha[[1]]-1)*log(continuousChoice[[1]])
    for(j in 2:length(V)){
      if(!anyNA(minConsumption)){
        tmp <- continuousChoice[[j]] - (continuousChoice[[j]] >= minConsumption[[j]])*minConsumption[[j]]
        V[[j]] = V[[j]] + avail[[j]]*((alpha[[j]]-1)*log((tmp/gamma[[j]]) + 1) - log(cost[[j]]))
      } else {
        V[[j]] = V[[j]] + avail[[j]]*((alpha[[j]]-1)*log((continuousChoice[[j]]/gamma[[j]]) + 1) - log(cost[[j]]+(1-avail[[j]])))
      }
    }
    
    
    #PART 1: Jacobian determinant
    fi = mapply(function(a, mc, g) (1-a)/(mc+g), alpha[-1], continuousChoice[-1], gamma[-1], SIMPLIFY=FALSE)
    term1_1 = mapply(function(f,d) f^d, fi, discrete_choice[-1], SIMPLIFY=FALSE)
    term1_1 = exp(Reduce("+", lapply(term1_1, log)))
    term1_1 = term1_1 * (1 - alpha[[1]]) / continuousChoice[[1]]
    term1_2 = mapply(function(c,f,d) c/f*d, cost[-1], fi, discrete_choice[-1], SIMPLIFY=FALSE)
    term1_2 = Reduce("+", term1_2) + continuousChoice[[1]]/(1 - alpha[[1]])
    term1   = term1_1*term1_2

    # Nest denominators
    nestDenom <- list()
    for(s in 1:nNests){
      alts <- which(as.vector(omega[s,])>0)
      nestDenom[[s]] <- Reduce("+", mapply(function(v,t) exp(v/t), V[alts], mdcnevNestsAlt[alts], SIMPLIFY=FALSE))
    }

    #PART 2
    term2_numerator = mapply(function(v, t, d) d*(v/t), V, mdcnevNestsAlt, discrete_choice, SIMPLIFY=FALSE)
    term2_numerator = exp(Reduce("+", term2_numerator))
    term2_denom = Reduce("*", mapply(function(nd, qs) nd^qs, nestDenom, q, SIMPLIFY=FALSE)) # if no product for a nest is consumed, then qs=0 and it is excluded
    term2 = term2_numerator / term2_denom

    #Term 3, part inside square brackets ==> this is done at the level of the whole sample:
    #Numerator
    term3_num   <- mapply(function(nd, t) nd^t, nestDenom, mdcnevNests, SIMPLIFY=FALSE)
    term3_denom <- Reduce("+", term3_num)
    term3SqBrac <- lapply(term3_num, "/", term3_denom) # David says: this wasn't here in the original code.
    tmp  <- term3SqBrac[[1]]
    pVec <- is.vector(tmp)
    pMat <- is.matrix(tmp)
    pCub <- (is.array(tmp) && length(dim(tmp))==3)
    rm(tmp)

    #TERM4
    r_current_combo=vector("double", nNests)

    
    
    #NOW WE LOOK AT THE SUMS OUTSIDE THE CURLY BRACKETS

    term345_total=vector("list", nObs)

    n=1
    while(n<(nObs+1)) {
      q_person <- do.call(c,lapply(q, function(qq) qq[n]))
      # create index of combinations for sums
      x = vector("list", sum(q_person>0))
      chosen_nests= vector("double", sum(q_person>0))

      i <- 1
      for(s in 1:nNests) if(q_person[s]>0){
        x[[i]] = 1:q_person[s]
        chosen_nests[i] = s
        i = i + 1
      }; rm(i, s)

      # create combinations
      sum_combo = expand.grid(x)

      # in some of the sums or products, there are only M different r vectors, but we use K different ones as we always go over all K nests, and just exclude those where q=0

      # create a new vector to contain current combination
      #r_current_combo=vector("double",S)
      # contribution to term345 is now calculated for every combination

      term345 = vector("list", nrow(sum_combo))

      for(kk in 1:nrow(sum_combo)){
        i <- 1
        for(s in 1:nNests){
          if(q_person[s]!=0){
            r_current_combo[s] = sum_combo[kk,i]
            i = i + 1
          } else r_current_combo[s] <- 0
        }
        term345[[kk]] = 1 #initialise to 1 so we can start multiplying
        for(s in 1:nNests) if(q_person[s]>0){
          if(pVec) tmp <- term3SqBrac[[s]][n]
          if(pMat) tmp <- term3SqBrac[[s]][n,]
          if(pCub) tmp <- term3SqBrac[[s]][n,,,drop=FALSE]
          term345[[kk]] = term345[[kk]]*tmp^(q_person[s] - r_current_combo[s] + 1)
        }
        term345[[kk]] = term345[[kk]]*term4final(mdcnevNests, r_current_combo, q_person)
        term345[[kk]] = term345[[kk]]*term5(r_current_combo, q_person)
      }
      term345_total[[n]] = Reduce("+", term345)
      n=n+1
    }

    if(pVec) term345_total <- do.call(c, term345_total)
    if(pMat) term345_total <- do.call(rbind, term345_total)
    if(pCub){
      tmp <- array(0, dim=c(nObs, dim(term345_total[[1]])[2:3]))
      for(n in 1:nObs) tmp[n,,] <- term345_total[[n]]
    }

    P <- term1*term2*term345_total
    if(pVec) P[!rows]   <- 1
    if(pMat) P[!rows,]  <- 1
    if(pCub) P[!rows,,] <- 1

    # make the chosen unavailable alternatives have a likelihood of zero
    choseUnavail <- mapply(function(m,a) m & !a, discrete_choice, avail, SIMPLIFY=TRUE)
    choseUnavail <- rowSums(choseUnavail)>0
    if(is.vector(P)) P[choseUnavail]   <- 0
    if(is.matrix(P)) P[choseUnavail,]  <- 0
    if(is.array(P) && length(dim(P))==3) P[choseUnavail,,] <- 0
    
    return( P )
  }


  # ############################## #
  #### functionality="output" ####
  # ############################## #

  if(functionality=="output"){
    # Calculate likelihood
    P <- apollo_mdcnev(mdcnev_settings,functionality="estimate")
    # Useful values
    omega <- mdcnevStructure
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    nNests<- length(mdcnevNests)
    avail_set <- FALSE
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      for(i in which(!rows)){
        continuousChoice[[j]][1] <- budget[i]
        for(j in 2:nAlts) continuousChoice[[j]][i] <- 0
      }
    }
    # Create availability if needed
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), alternatives)
    }
    discrete_choice=list()
    for(k in 1:nAlts) discrete_choice[[alternatives[k]]] = (continuousChoice[[k]]>0)*1

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
    colnames(choicematrix) = c(outsideName, names(V)[2:length(V)])
    colnames(mdcnevStructure) <- c(outsideName, names(V)[2:length(V)])
    rownames(mdcnevStructure) <- names(mdcnevNests)
    content <- list(round(choicematrix,2))
    if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                         "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    content[[length(content) + 1]] <- "Structure for MDCNEV model component:"
    content[[length(content) + 1]] <- mdcnevStructure
    apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    apollo_addLog("Overview of choices for MDCEV model component:", content, apolloLog)
    apollo_addLog(title   = paste0("Final structure for MDCNEV model component \"",componentName,"\":"), 
                  content = mdcnevStructure, apolloLog, 
                  book    = 2)
    apollo_reportModelTypeLog(modelType="MDCNEV", apolloLog)

    return(P)
  }

  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction"){
    # Store useful values
    nObs  <- length(continuousChoice[[1]])
    nAlts <- length(V)
    if(length(rows)==1 && rows=="all") rows <- rep(TRUE, nObs) else {
      if(length(budget)==length(continuousChoice[[1]])) continuousChoice[[1]][!rows] = budget[!rows] else {
        continuousChoice[[1]][!rows] = budget
      }
      for(j in 2:nAlts) continuousChoice[[j]][!rows] = 0
    }

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
      nAlts <- length(continuousChoice)
      nNests <- nrow(mdcnevStructure)

      corr_matrix=matrix(0,nAlts,nAlts)
      for(n in 1:nNests){
        alts <- which(mdcnevStructure[n,]==1)
        corr_matrix[alts, alts] <- 1 - mdcnevNests[[n]]^2
      }
      diag(corr_matrix) <- 1
      covMat <- corr_matrix  # covariance matrix = corr_matrix as all std dev are 1

      ND = 250
      Edraws = mnormt::rmnorm(n=ND*nObs, mean=rep(0, nAlts), varcov=covMat)
      Edraws = stats::pnorm(Edraws)
      Edraws = -log(-log(Edraws))


      k=1
      cat("\n0%")
      while(k<(ND+1))
      {
        draws <- sigma*Edraws[(1 + (k-1)*nObs):(k*nObs), ]

        V_trans=V

        V_trans[[1]] = exp(V_trans[[1]] + draws[,1])

        j=2
        while(j<=length(V)){
          V_trans[[j]] = exp(V_trans[[j]] + draws[,j])
          V_trans[[j]] = V_trans[[j]]/cost[[j]]
          V_trans[[j]] = V_trans[[j]]*avail[[j]]
          j=j+1
        }

        for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
        for(i in 1:length(alpha)) if(length(alpha[[i]])==1) alpha[[i]] <- rep(alpha[[i]], nObs)
        for(i in 1:length(gamma)) if(length(gamma[[i]])==1) gamma[[i]] <- rep(gamma[[i]], nObs)
        for(i in 1:length(cost)) if(length(cost[[i]])==1) cost[[i]] <- rep(cost[[i]], nObs)

        # need a check here to see if we're using a given column for example of the cube

        # turn into matrices
        V_trans_use = Reduce("cbind", V_trans)
        alpha_use   = Reduce("cbind", alpha)
        gamma_use   = Reduce("cbind", gamma)
        avail_use   = Reduce("cbind", avail)
        cost_use    = Reduce("cbind", cost)

        consumption_total[,,k] = forecasting_subroutine(V_trans_use,alpha_use,gamma_use,avail_use,cost_use,budget,continuousChoice)

        # increment draws
        if(k%%round(ND/10,0)==0) if(k==round(ND/2,0)) cat("50%") else if(k==ND) cat("100%") else cat(".")
        k=k+1
      }
      cat("\n")
      return(consumption_total)
    }

    forecasting_subroutine=function(V_trans, alpha, gamma, avail, cost, budget, continuousChoice)
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
          mdcnevNests_1=((budget[i]+sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use)))
          mdcnevNests_21=(V_trans[i,1])^(1/(1-alpha[i,1]))
          mdcnevNests_22=sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use*(V_trans[i,2:NALT])^(1/(1-alpha[i,1])))
          mdcnevNests=(mdcnevNests_1/(mdcnevNests_21+mdcnevNests_22))^(alpha[i,1]-1)
          if( sum(mdcnevNests<V_trans[i,2:NALT]) < M )
          {
            #step3
            consumption_outside_1=(V_trans[i,1])^(1/(1-alpha[i,1]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
            consumption_outside_2=mdcnevNests_21+mdcnevNests_22
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
              mdcnevNests_22=sum(cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]*use*(V_trans[i,2:NALT])^(1/(1-alpha[i,1])))
              consumption_outside_1=(V_trans[i,1])^(1/(1-alpha[i,1]))*(budget[i]+sum(use*cost[i,(2:ncol(cost))]*gamma[i,(2:ncol(gamma))]))
              consumption_outside_2=mdcnevNests_21+mdcnevNests_22
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

      dimnames(consumption)[[2]]=alternatives

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
