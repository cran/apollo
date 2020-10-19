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
#' @importFrom mnormt rmnorm
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @export
apollo_mdcnev <- function(mdcnev_settings,functionality){
  ### Set or extract componentName
  modelType = "MDCNEV"
  if(is.null(mdcnev_settings[["componentName"]])){
    mdcnev_settings[["componentName"]] = ifelse(!is.null(mdcnev_settings[['componentName2']]),
                                                mdcnev_settings[['componentName2']], modelType)
    test <- functionality=="validate" && mdcnev_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 mdcnev_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, mdcnev_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", mdcnev_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(mdcnev_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load mdcnev_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(mdcnev_settings$componentName, "_settings")]]
    # If there is no V inside the loaded mdcnev_settings, restore the one received as argument
    if(is.null(tmp$V    )) tmp$V     <- mdcnev_settings$V    
    if(is.null(tmp$alpha)) tmp$alpha <- mdcnev_settings$alpha
    if(is.null(tmp$gamma)) tmp$gamma <- mdcnev_settings$gamma
    if(is.null(tmp$sigma)) tmp$sigma <- mdcnev_settings$sigma
    if(is.null(tmp[["mdcnevNests"]])) tmp$mdcnevNests <- mdcnev_settings$mdcnevNests
    mdcnev_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    mdcnev_settings <- apollo_preprocess(inputs = mdcnev_settings, modelType, 
                                         functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print("No C++ optimisation available for MDCNEV")
    mdcnev_settings$probs_MDCNEV <- function(mdcnev_settings){
      # Set utility of unavailable alternatives and excluded rows to 0 to avoid numerical issues
      mdcnev_settings$V <- mapply(function(v,a) apollo_setRows(v, !a , 0), 
                                  mdcnev_settings$V, mdcnev_settings$avail, SIMPLIFY=FALSE)
      # Compute V
      mdcnev_settings$V[[1]] = (mdcnev_settings$alpha[[1]]-1)*log(mdcnev_settings$continuousChoice[[1]])
      for(j in 2:mdcnev_settings$nAlt){
        if(!anyNA(mdcnev_settings$minConsumption)){
          tmp <- mdcnev_settings$continuousChoice[[j]] - 
            (mdcnev_settings$continuousChoice[[j]] >= mdcnev_settings$minConsumption[[j]])*mdcnev_settings$minConsumption[[j]]
          mdcnev_settings$V[[j]] = mdcnev_settings$V[[j]] + 
            mdcnev_settings$avail[[j]]*((mdcnev_settings$alpha[[j]]-1)*log(tmp/mdcnev_settings$gamma[[j]] + 1) - 
                                          mdcnev_settings$log(mdcnev_settings$cost[[j]]))
        } else {
          mdcnev_settings$V[[j]] = mdcnev_settings$V[[j]] + 
            mdcnev_settings$avail[[j]]*((mdcnev_settings$alpha[[j]]-1)*log(mdcnev_settings$continuousChoice[[j]]/mdcnev_settings$gamma[[j]] + 1) - 
                                          log(mdcnev_settings$cost[[j]]+(1-mdcnev_settings$avail[[j]])))
        }
      }
      #PART 1: Jacobian determinant
      fi = mapply(function(a, mc, g) (1-a)/(mc+g), 
                  mdcnev_settings$alpha[-1], mdcnev_settings$continuousChoice[-1], 
                  mdcnev_settings$gamma[-1], SIMPLIFY=FALSE)
      term1_1 = mapply(function(f,d) f^d, fi, mdcnev_settings$discrete_choice[-1], SIMPLIFY=FALSE)
      term1_1 = exp(Reduce("+", lapply(term1_1, log)))
      term1_1 = term1_1 * (1 - mdcnev_settings$alpha[[1]]) / mdcnev_settings$continuousChoice[[1]]
      term1_2 = mapply(function(c,f,d) c/f*d, 
                       mdcnev_settings$cost[-1], fi, mdcnev_settings$discrete_choice[-1], SIMPLIFY=FALSE)
      term1_2 = Reduce("+", term1_2) + mdcnev_settings$continuousChoice[[1]]/(1 - mdcnev_settings$alpha[[1]])
      term1   = term1_1*term1_2
      # Nest denominators
      altsInNest = list() # called theta in old code
      for(k in 1:mdcnev_settings$nAlt) altsInNest[[k]] <- 
        mdcnev_settings$mdcnevNests[[which(as.vector(mdcnev_settings$mdcnevStructure[,k])>0)]]
      nestDenom <- list()
      for(s in 1:mdcnev_settings$nNests){
        alts <- which(as.vector(mdcnev_settings$mdcnevStructure[s,])>0)
        nestDenom[[s]] <- Reduce("+", mapply(function(v,t) exp(v/t), 
                                             mdcnev_settings$V[alts], altsInNest[alts], SIMPLIFY=FALSE))
      }
      #PART 2
      term2_numerator = mapply(function(v, t, d) d*(v/t), mdcnev_settings$V, 
                               altsInNest, mdcnev_settings$discrete_choice, SIMPLIFY=FALSE)
      term2_numerator = exp(Reduce("+", term2_numerator))
      term2_denom = Reduce("*", mapply(function(nd, qs) nd^qs, nestDenom, mdcnev_settings$q, SIMPLIFY=FALSE)) # if no product for a nest is consumed, then qs=0 and it is excluded
      term2 = term2_numerator / term2_denom
      #Term 3, part inside square brackets ==> this is done at the level of the whole sample:
      #Numerator
      term3_num   <- mapply(function(nd, t) nd^t, nestDenom, mdcnev_settings$mdcnevNests, SIMPLIFY=FALSE)
      term3_denom <- Reduce("+", term3_num)
      term3SqBrac <- lapply(term3_num, "/", term3_denom) # David says: this wasn't here in the original code.
      tmp  <- term3SqBrac[[1]]
      pVec <- is.vector(tmp)
      pMat <- is.matrix(tmp)
      pCub <- (is.array(tmp) && length(dim(tmp))==3)
      rm(tmp)
      #TERM4
      r_current_combo=vector("double", mdcnev_settings$nNests)
      #NOW WE LOOK AT THE SUMS OUTSIDE THE CURLY BRACKETS
      term345_total=vector("list", mdcnev_settings$nObs)
      for(n in 1:mdcnev_settings$nObs) {
        q_person <- do.call(c,lapply(mdcnev_settings$q, function(qq) qq[n]))
        # create index of combinations for sums
        x = vector("list", sum(q_person>0))
        chosen_nests= vector("double", sum(q_person>0))
        i <- 1
        for(s in 1:mdcnev_settings$nNests) if(q_person[s]>0){
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
          for(s in 1:mdcnev_settings$nNests){
            if(q_person[s]!=0){
              r_current_combo[s] = sum_combo[kk,i]
              i = i + 1
            } else r_current_combo[s] <- 0
          }
          term345[[kk]] = 1 #initialise to 1 so we can start multiplying
          for(s in 1:mdcnev_settings$nNests) if(q_person[s]>0){
            if(pVec) tmp <- term3SqBrac[[s]][n]
            if(pMat) tmp <- term3SqBrac[[s]][n,]
            if(pCub) tmp <- term3SqBrac[[s]][n,,,drop=FALSE]
            term345[[kk]] = term345[[kk]]*tmp^(q_person[s] - r_current_combo[s] + 1)
          }
          term345[[kk]] = term345[[kk]]*mdcnev_settings$term4final(mdcnev_settings$mdcnevNests, r_current_combo, q_person)
          term345[[kk]] = term345[[kk]]*mdcnev_settings$term5(r_current_combo, q_person)
        }
        term345_total[[n]] = Reduce("+", term345)
      }
      if(pVec) term345_total <- do.call(c, term345_total)
      if(pMat) term345_total <- do.call(rbind, term345_total)
      if(pCub){
        tmp <- array(0, dim=c(mdcnev_settings$nObs, dim(term345_total[[1]])[2:3]))
        for(n in 1:mdcnev_settings$nObs) tmp[n,,] <- term345_total[[n]]
      }
      P <- term1*term2*term345_total
      # make the chosen unavailable alternatives have a likelihood of zero
      choseUnavail <- mapply(function(m,a) m & !a, mdcnev_settings$discrete_choice, mdcnev_settings$avail, SIMPLIFY=TRUE)
      choseUnavail <- rowSums(choseUnavail)>0
      if(is.vector(P)) P[choseUnavail]   <- 0
      if(is.matrix(P)) P[choseUnavail,]  <- 0
      if(is.array(P) && length(dim(P))==3) P[choseUnavail,,] <- 0
      return( P )
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    apollo_beta <- tryCatch(get("apollo_beta", envir=parent.frame(), inherits=TRUE),
                            error=function(e) return(NULL))
    test <- !is.null(apollo_beta) && (functionality %in% c("preprocess", "gradient"))
    test <- test && all(sapply(mdcnev_settings$V, is.function))
    test <- test && is.function(mdcnev_settings$alpha)
    test <- test && is.function(mdcnev_settings$gamma)
    test <- test && is.function(mdcnev_settings$sigma)
    test <- test && apollo_inputs$apollo_control$analyticGrad
    mdcnev_settings$gradient <- FALSE
    if(test){
      mdcnev_settings$dV     <- apollo_dVdB(apollo_beta, apollo_inputs, mdcnev_settings$V)
      mdcnev_settings$dAlpha <- apollo_dVdB(apollo_beta, apollo_inputs, mdcnev_settings$alpha)
      mdcnev_settings$dGamma <- apollo_dVdB(apollo_beta, apollo_inputs, mdcnev_settings$gamma)
      mdcnev_settings$dSigma <- apollo_dVdB(apollo_beta, apollo_inputs, list(dSigma=mdcnev_settings$sigma))[[1]]
      #mdcnev_settings$gradient <- !is.null(mdcnev_settings$dV)
    }; rm(test)
    
    # Return mdcnev_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      mdcnev_settings$V     <- NULL
      mdcnev_settings$alpha <- NULL
      mdcnev_settings$gamma <- NULL
      mdcnev_settings$sigma <- NULL
      mdcnev_settings$mdcnevNests <- NULL
      return(mdcnev_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute V, alpha, gamma and sigma (makes sure we are now working with vectors/matrices/arrays and not functions)
  testV <- any(sapply(mdcnev_settings$V    , is.function))
  testA <- any(sapply(mdcnev_settings$alpha, is.function))
  testG <- any(sapply(mdcnev_settings$gamma, is.function))
  testS <- is.function(mdcnev_settings$sigma)
  testN <- any(sapply(mdcnev_settings$mdcnevNests, is.function))
  if(testV) mdcnev_settings$V     = lapply(mdcnev_settings$V    , function(f) if(is.function(f)) f() else f )
  if(testA) mdcnev_settings$alpha = lapply(mdcnev_settings$alpha, function(f) if(is.function(f)) f() else f )
  if(testG) mdcnev_settings$gamma = lapply(mdcnev_settings$gamma, function(f) if(is.function(f)) f() else f )
  if(testS) mdcnev_settings$sigma = mdcnev_settings$sigma()
  if(testN) mdcnev_settings$mdcnevNests = lapply(mdcnev_settings$mdcnevNests, function(f) if(is.function(f)) f() else f )
  rm(testV, testA, testG, testS, testN)
  mdcnev_settings$V     <- lapply(mdcnev_settings$V    , function(v) if(is.matrix(v) && ncol(v)==1) as.vector(v) else v)
  mdcnev_settings$alpha <- lapply(mdcnev_settings$alpha, function(a) if(is.matrix(a) && ncol(a)==1) as.vector(a) else a)
  mdcnev_settings$gamma <- lapply(mdcnev_settings$gamma, function(g) if(is.matrix(g) && ncol(g)==1) as.vector(g) else g)
  mdcnev_settings$mdcnevNests <- lapply(mdcnev_settings$mdcnevNests, function(g) if(is.matrix(g) && ncol(g)==1) as.vector(g) else g)
  if(is.matrix(mdcnev_settings$sigma) && ncol(mdcnev_settings$sigma)==1) mdcnev_settings$sigma <- as.vector(mdcnev_settings$sigma)
  
  ### Deal with outside
  if(mdcnev_settings$hasOutside){
    # Add gamma outside, if missing
    if(is.null(mdcnev_settings$gamma$outside)) mdcnev_settings$gamma$outside <- 1
    # Replace mdcnev_settings$outside by "outside" in V, alpha, and gamma
    tmp <- which(names(mdcnev_settings$V    )==mdcnev_settings$output); if(length(tmp)>0) names(mdcnev_settings$V    )[tmp] <- "outside"
    tmp <- which(names(mdcnev_settings$alpha)==mdcnev_settings$output); if(length(tmp)>0) names(mdcnev_settings$alpha)[tmp] <- "outside"
    tmp <- which(names(mdcnev_settings$gamma)==mdcnev_settings$output); if(length(tmp)>0) names(mdcnev_settings$gamma)[tmp] <- "outside"
    rm(tmp)
  }
  ### Reorder V, alpha and gamma if necessary
  if( any(mdcnev_settings$alternatives != names(mdcnev_settings$V    )) ) mdcnev_settings$V     <- mdcnev_settings$V[mdcnev_settings$alternatives]
  if( any(mdcnev_settings$alternatives != names(mdcnev_settings$alpha)) ) mdcnev_settings$alpha <- mdcnev_settings$alpha[mdcnev_settings$alternatives]
  if( any(mdcnev_settings$alternatives != names(mdcnev_settings$gamma)) ) mdcnev_settings$gamma <- mdcnev_settings$gamma[mdcnev_settings$alternatives]
  
  ### Drop rows if neccesary
  if(!all(mdcnev_settings$rows)){
    mdcnev_settings$V     <- lapply(mdcnev_settings$V    , apollo_keepRows, r=mdcnev_settings$rows)
    mdcnev_settings$alpha <- lapply(mdcnev_settings$alpha, apollo_keepRows, r=mdcnev_settings$rows)
    mdcnev_settings$gamma <- lapply(mdcnev_settings$gamma, apollo_keepRows, r=mdcnev_settings$rows)
    mdcnev_settings$sigma <- apollo_keepRows(mdcnev_settings$sigma, r=mdcnev_settings$rows)
    mdcnev_settings$mdcnevNests <- lapply(mdcnev_settings$mdcnevNests, apollo_keepRows, r=mdcnev_settings$rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #

  if(functionality=="validate"){
    
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(mdcnev_settings, modelType, 
                                                                   functionality, apollo_inputs)

    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(mdcnev_settings, modelType, apollo_inputs)
    
    testL = mdcnev_settings$probs_MDCNEV(mdcnev_settings)
    if(any(!mdcnev_settings$rows)) testL <- apollo_insertRows(testL, mdcnev_settings$rows, 1)
    if(all(testL==0)) stop('All observations have zero probability at starting value for model component "', 
                           mdcnev_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', 
                                          mdcnev_settings$componentName,'"'))
    return(invisible(testL))
  }
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #

  if(functionality=="zero_LL"){
    P <- rep(NA, mdcnev_settings$nObs)
    P[!mdcnev_settings$rows] <- 1
    return(P)
  }

  # ################################################## #
  #### functionality="estimate/conditionals/output" ####
  # ################################################## #
  if(functionality %in% c("estimate", "conditionals", "output")){
    L <- mdcnev_settings$probs_MDCNEV(mdcnev_settings)
    if(any(!mdcnev_settings$rows)) L <- apollo_insertRows(L, mdcnev_settings$rows, 1)
    return(L)
  }

  # ################################ #
  #### functionality="prediction" ####
  # ################################ #

  if(functionality=="prediction"){
    # Change name to mdcnev_settings to "s"
    s <- mdcnev_settings
    rm(mdcnev_settings)
    
    # Check that sigma is not random
    if(!is.vector(s$sigma)) stop('Forecasting not available for random sigma in models component ', 
                                 s$componentName)
    
    # Generate draws for correlated gumbel error components
    set.seed(99)
    corr_matrix = matrix(0, s$nAlt, s$nAlt)
    for(n in 1:s$nNests){
      alts <- which(s$mdcnevStructure[n,]==1)
      corr_matrix[alts, alts] <- 1 - s$mdcnevNests[[n]]^2
    }; rm(alts)
    diag(corr_matrix) <- 1 # covariance matrix = corr_matrix as all std dev are 1
    tmp1 <- mnormt::rmnorm(n=s$nRep*s$nObs, mean=rep(0, s$nAlt), varcov=corr_matrix)
    tmp1 <- stats::pnorm(tmp1)
    tmp1 = -log(-log(tmp1))
    epsL <- vector(mode="list", length=s$nRep)
    tmp2 <- (0:(s$nObs-1))*s$nRep
    for(r in 1:s$nRep) epsL[[r]] <- tmp1[r+tmp2,]#split(tmp1[r+tmp2,], f=rep(1:s$nAlt, each=s$nObs))
    rm(tmp1, tmp2)
    
    # Initialise output matrices
    Xm <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    Xv <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    Mm <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    Mv <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    X  <- matrix(0, nrow=s$nObs, ncol=s$nAlt)
    
    # Tools to deal with draws
    extractDraw <- function(b, iInter, iIntra){
      if(!is.list(b)){
        if(is.vector(b) && length(b)==1) return(rep(b, s$nObs))
        if(is.vector(b)) return(b)
        if(is.matrix(b)) return(b[,iInter])
        if(is.array(b) && length(dim(b))==3) return(b[,iInter,iIntra])
      } else {
        ans <- lapply(b, extractDraw, iInter=iInter, iIntra=iIntra)
        ans <- do.call(cbind, ans)
        return(ans)
      }
    }
    if(anyNA(apollo_inputs$apollo_draws)){
      nInter <- 0; nIntra <- 0
    } else {
      nInter <- apollo_inputs$apollo_draws$interNDraws
      nIntra <- apollo_inputs$apollo_draws$intraNDraws
    }
    if(nInter==0) nInter <- 1
    if(nIntra==0) nIntra <- 1
    step1 <- ceiling(s$nRep/10)
    step2 <- ceiling(s$nRep/2)
    
    # Expand avail and cost
    avail <- sapply(s$avail, function(a) if(length(a)==1) rep(a, s$nObs) else a)
    cost  <- sapply(s$cost , function(x) if(length(x)==1) rep(x, s$nObs) else x)
    
    ### Simulate prediction
    if(!is.null(apollo_inputs$silent)) silent <- apollo_inputs$silent else silent <- FALSE
    if(!silent) cat("0%")
    for(r in 1:s$nRep){
      X    <- 0*X
      for(iInter in 1:nInter){
        Xintra <- 0*X
        for(iIntra in 1:nIntra){
          # Extract appropiate values
          phiP <- extractDraw(s$sigma, iInter, iIntra)*epsL[[r]]
          phiP <- exp(extractDraw(s$V, iInter, iIntra) + phiP)*avail/cost
          gamma<- extractDraw(s$gamma, iInter, iIntra)
          alpha<- extractDraw(s$alpha, iInter, iIntra)
          # Calculate prediction
          for(i in 1:s$nObs){
            p  <- cost[i,2:s$nAlt]
            b  <- s$budget[i]
            g  <- gamma[i,2:s$nAlt]
            a0 <- alpha[i,1]
            ak <- alpha[i,2:s$nAlt]
            ph0<- phiP[i,1] 
            phk<- phiP[i,2:s$nAlt] 
            orderofV = rank(-phk)
            M = 1
            stopping = FALSE
            while(!stopping){ # nAlt
              use = orderofV<M
              #step2
              lambda_1  = b + sum(p*g*use)
              lambda_21 = ph0^(1/(1-a0))
              lambda_22 = sum(p*g*use*phk^(1/(1-a0)))
              lambda_2  = lambda_21 + lambda_22
              lambda    = (lambda_1/lambda_2)^(a0-1)
              if( M > sum(phk>lambda) || M > sum(phk>0)){
                #step3
                x0_1 = lambda_21*lambda_1
                Xintra[i,1] = Xintra[i,1] + x0_1/lambda_2 # maybe X = X + ... to deal with draws
                xk_1 = phk^(1/(1-ak))*lambda_1
                Xintra[i,2:s$nAlt] = Xintra[i,2:s$nAlt] + use*(xk_1/lambda_2 - 1)*g  # maybe X = X + ... to deal with draws
                stopping = TRUE
              } else M <- M + 1 # step 4
            } # end of "while"
          } # end of "for" loop over observations
        } # end of "for" loop over intra draws
        X <- X + Xintra/nIntra # CHECK IF THIS IS CORRECT!!!
      } # end of "for" loop over inter draws
      X <- X/nInter # CHECK IF THIS IS CORRECT!!!
      # Store prediction (maybe this needs to be nested in the for loops in a more clever way to address for draws)
      Xm <- Xm + X/s$nRep
      Mm <- Mm + (X>0)/s$nRep
      Xv <- Xv + apply(X  , MARGIN=2, function(v) (v-mean(v))^2)/s$nRep
      Mv <- Mv + apply(X>0, MARGIN=2, function(m) (m-mean(m))^2)/s$nRep
      if(!silent){ if(r%%step2==0) cat(round(100*r/s$nRep,0), "%", sep="") else { if(r%%step1==0) cat(".") } }
    } # end of "for" loop over repetitions
    if(!silent) cat("\n")
    
    ### Prepare output
    out <- cbind(Xm, sqrt(Xv), Mm, sqrt(Mv))
    out <- apollo_insertRows(out, s$rows, NA)
    #colnames(out) <- paste( names(s$continuousChoice), rep(c("ECont", "sdCont", "EDisc", "sdDisc"), each=s$nAlt), sep="_")
    colnames(out) <- paste( names(s$continuousChoice), rep(c("cont_mean", "cont_sd", "disc_mean", "disc_sd"), each=s$nAlt), sep="_")
    return(out)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(mdcnev_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(mdcnev_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
}
