#' MDC model with exogenous budget
#' 
#' Calculates the likelihood function of the MDC model with exogenous budget. Can also predict and validate inputs.
#' 
#' This model extends the traditional multiple discrete-continuous (MDC) framework by (i) making the 
#' marginal utility of the outside good deterministic, and (ii) including complementarity and 
#' substitution in the model formulation. See the following working paper for more details:
#' 
#' Palma, D. & Hess, S. (2022) Extending the Multiple Discrete Continuous (MDC) modelling 
#' framework to consider complementarity, substitution, and an unobserved budget. Transportation 
#' Reserarch 161B, 13 - 35. https://doi.org/10.1016/j.trb.2022.04.005
#' 
#' @param emdc_settings List of settings for the model. It includes the following.
#'                        \itemize{
#'                          \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount consumed of each inside good. Outside good must not be included. Can also be called "X".
#'                          \item \strong{\code{budget}}: Numeric vector. Budget. Must be bigger that the expenditure on all inside goods. Can also be called "B".
#'                          \item \strong{\code{avail}}: Named list of numeric vectors. Availability of each product. Can also be called "A".
#'                          \item \strong{\code{utilityOutside}}: Numeric vector (or matrix or array). Shadow price of the budget. Must be normalised to 0 for at least one individual. Default is 0 for every observation. Can also be called "V0".
#'                          \item \strong{\code{utilities}}: Named list of numeric vectors (or matrices or arrays). Base utility of each product. Can also be called "V".
#'                          \item \strong{\code{gamma}}: Named list of numeric vectors. Satiation parameter of each product.
#'                          \item \strong{\code{delta}}: Lower triangular numeric matrix, or list of lists. Complementarity/substitution parameter.
#'                          \item \strong{\code{cost}}: Named list of numeric vectors. Price of each product.
#'                          \item \strong{\code{sigma}}: Numeric vector or scalar. Standard deviation of the error term. Default is one.
#'                          \item \strong{\code{nRep}}: Scalar positive integer. Number of repetitions used when prediction
#'                          \item \strong{\code{tol}}: Positive scalar. Tolerance of the prediction algorithm.
#'                          \item \strong{\code{timeLimit}}: Positive scalar. Maximum amount of seconds the optimiser can spend calculating a prediction before setting it to NA.
#'                        }
#' @param functionality Character. Either "validate", "zero_LL", "estimate", "conditionals", "raw", "output" or "prediction"
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
#' @export
#' @importFrom stats qnorm
#' @importFrom Rsolnp solnp
apollo_emdc1 <- function(emdc_settings, functionality="estimate"){
  # Rename input if necessary
  map <- c(X = "continuousChoice", B = "budget", A = "avail", 
           V0= "utilityOutside",   V = "utilities")
  for(i in 1:length(map)) if(!is.null(emdc_settings[[map[i]]])){
    emdc_settings[[names(map)[i]]] <- emdc_settings[[map[i]]]
    emdc_settings[[map[i]]]        <- NULL
  }; rm(i, map)
  
  # Check input
  mandatory <- c("X", "B", "A", "V", "gamma", "delta", "cost")
  optional  <- list(V0=0, nRep=50, tol=0.1, sigma=1, timeLimit=5*60)
  test <- mandatory %in% names(emdc_settings)
  if(!all(test)) stop('Mandatory setting(s) "', paste0(mandatory[!test], collapse='", "'),
                      '" are missing from "emdc_settings".')
  for(i in names(optional)) if(!(i %in% names(emdc_settings))) emdc_settings[[i]] <- optional[[i]]
  test <- is.list(emdc_settings$X) && all(sapply(emdc_settings$X, is.vector)) && !is.null(names(emdc_settings$X))
  if(!test) stop("'X' should be a named list of numeric vectors.")
  nAlt <- length(emdc_settings$X)
  nObs <- max(sapply(emdc_settings$X, length))
  test <- is.vector(emdc_settings$B) && is.numeric(emdc_settings$B)
  test <- test && (length(emdc_settings$B) %in% c(1,nObs)) && all(emdc_settings$B>0)
  if(!test) stop('Setting "B" must be a numeric vector of non-negative values with 1 or nObs elements.')
  test <- is.list(emdc_settings$A) && all(sapply(emdc_settings$A, is.vector))
  test <-  test && all(sapply(emdc_settings$A, length) %in% c(1, nObs))
  if(!test) stop('Setting "A" should be a list of numeric vectors, each of length 1 or nObs.')
  test <- is.vector(emdc_settings$V0) | is.array(emdc_settings$V0)
  test <- test && ifelse(is.array(emdc_settings$V0), dim(emdc_settings$V0), length(emdc_settings$V0)) %in% c(1,nObs)
  if(!test) stop("'V0' should be a numeric vector, matrix or array, with 1 or nObs rows.")
  test <- is.list(emdc_settings$V) && all(sapply(emdc_settings$V, function(v) is.vector(v) | is.array(v)))
  if(!test) stop("'V' should be a list of numeric vectors, matrices or arrays.")
  test <- is.list(emdc_settings$gamma) && all(sapply(emdc_settings$gamma, function(g) is.vector(g) | is.array(g)))
  if(!test) stop("'gamma' should be a list of numeric vectors, matrices or arrays.")
  test <- is.list(emdc_settings$cost) && all(sapply(emdc_settings$cost, function(g) is.vector(g) | is.array(g)))
  test <- all(sapply(emdc_settings$cost, is.numeric))
  if(!test) stop("'cost' should be a list of numeric vectors.")
  test <- is.vector(emdc_settings$nRep) && length(emdc_settings$nRep)==1 && emdc_settings$nRep>0
  if(!test) stop("'nRep' should be a positive scalar integer.")
  test <- is.vector(emdc_settings$tol) && length(emdc_settings$tol)==1 && emdc_settings$tol>0
  if(!test) stop("'tol' should be a positive numeric scalar.")
  test <- is.vector(emdc_settings$sigma) && is.numeric(emdc_settings$sigma) && length(emdc_settings$sigma) %in% c(1,nObs)
  if(!test) stop('Argument "sigma" must be a numerical scalar or a vector.')
  rm(mandatory, optional, test, i)
  
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE))
  
  # Copy variables from list to environment
  for(i in 1:length(emdc_settings)) assign(names(emdc_settings)[i], emdc_settings[[i]])
  
  # Create avail if necessary
  avail_set = FALSE
  if(length(A)==1 && A==1){
    A <- as.list(setNames(rep(1, length(X)), names(X)))
    avail_set <- TRUE
  }
  
  # Re-order arguments
  B     <- B # this prevents a NOTE when checking the package.
  V0    <- V0# this prevents a NOTE when checking the package.
  nRep  <- nRep #this prevents a NOTE when checking the package.
  sigma <- sigma# this prevents a NOTE when checking the package.
  timeLimit <- timeLimit # this prevents a NOTE when checking the package.
  A     <- A[names(X)]
  V     <- V[names(X)]
  gamma <- gamma[names(X)]
  cost  <- cost[names(X)]
  
  # Expand delta (lower triang matrix --> list of lists)
  # For example, for 3 products:
  # list( list(   0, d21, d31 ),
  #       list( d21,   0, d32 ),
  #       list( d31, d32,   0 ))
  test <- is.matrix(delta) && dim(delta)[1]==dim(delta)[2] && sum(delta[upper.tri(delta, diag=TRUE)])==0
  if(test){
    diag(delta) <- 0 # Important: makes it unnecessary to exclude i==j case in sums
    d    <- vector("list", nAlt)
    for(i in 1:nAlt){
      d[[i]] <- vector("list", nAlt)
      for(j in 1:nAlt) d[[i]][[j]] <- delta[max(i,j),min(i,j)]
    }
    delta <- d
    rm(d, i, j)
  } else {
    test <- is.list(delta) && all(sapply(delta, is.list)) && all(sapply(delta, length)==length(delta))
    if(test) for(i in 1:nAlt) test <- test && delta[[i]][[i]]==0
    if(test) for(i in 1:(nAlt-1)) for(j in (i+1):nAlt) test <- test && delta[[i]][[j]]==0
    if(test) for(i in 1:(nAlt-1)) for(j in (i+1):nAlt) delta[[i]][[j]] <- delta[[j]][[i]]
    #if(test) for(i in 2:nAlt) for(j in 1:(i-1)) test <- test && delta[[i]][[j]] = delta[[j]][[i]]
    if(!test) stop('"delta" should be a lower triangular matrix with zeros in the diagonal, or list of lists.')
  }
  
  
  # -------------- #
  #### VALIDATE ####
  # -------------- #
  if(functionality %in% c("validate")){
    for(a in names(X)){
      if(length(A[[a]])==1) A[[a]] <- rep(A[[a]], length(X[[a]]))
      if(any(A[[a]][X[[a]]>0]==0)) stop(paste0("Alternative ", a, " is chosen despite not being available"))
    }
    if(!(all(names(X)==names(A)) & all(names(A)==names(V)) & 
         all(names(V)==names(gamma)) & all(names(gamma)==names(cost)))) stop("Alternatives names are not the same across arguments")
    if(min(sapply(X, length))==0) stop("At least one element in X has length 0")
    if('outside' %in% names(X)) stop('The outside good must NOT be included in setting "X".')
    
    expenditure <- Reduce('+', mapply('*', X, cost, SIMPLIFY=FALSE))
    if(any(expenditure>=B)) stop('The budget must be bigger than the expenditure for all observations.')
    
    ### Create report
    r <- c("Times available", "Times chosen", "% chosen when avail.", 
           "Avg. consump. when avail.", "Avg. consump. when chosen")
    M <- matrix(0, nrow=5, ncol=length(V), dimnames=list(r, names(X)))
    for(a in names(X)){
      # Expand avail
      if(length(A[[a]])==1) A[[a]] <- rep(A[[a]], length(X[[a]]))
      M[1, a] <- sum(A[[a]])            # Times available
      M[2, a] <- sum(X[[a]]>0)          # Times chosen
      M[3, a] <- M[2, a]/M[1, a]*100    # % chosen when avail
      M[4, a] <- mean(X[[a]][A[[a]]>0]) # Avg. consump. when avail.
      M[5, a] <- mean(X[[a]][X[[a]]>0]) # Avg. consump. when chosen.
    }
    content <- list(round(M,2))
    if(any(M[5,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(M[3,]==100)) content[[length(content)+1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    
    testL <- apollo_emdc1(emdc_settings, functionality="estimate")
    if(all(testL==0)) stop("\nAll observations have zero probability at starting value for eMDCEV model component.")
    if(any(testL==0)) cat("\nSome observations have zero probability at starting value for eMDCEV model component.")
    return(invisible(testL))
  }
  
  # ------------- #
  #### ZERO LL ####
  # ------------- #
  if(functionality=="zero_LL"){
    ans <- rep(NA, length(X[[1]]))
    return(ans)
  }
  
  # ------------------------------------ #
  #### ESTIMATE, CONDITIONALS AND RAW ####
  # ------------------------------------ #
  if(functionality %in% c("estimate", "conditionals", "raw")){
    # Calculate outside good consumption
    expenditure <- Reduce('+', mapply('*', X, cost, SIMPLIFY=FALSE))
    x0 <- B - expenditure
    
    # Create useful variables
    M  <- lapply(X, function(x) x>0)
    Mm <- simplify2array(M)
    Mv <- rowSums(Mm)
    
    # Calculate C, D, E, G
    C <- mapply(function(xi, gi) 1/(xi + gi), X, gamma, SIMPLIFY=FALSE)
    D <- lapply(X, function(xi) exp(-xi))
    E <- vector(mode="list", nAlt)
    for(i in 1:nAlt) E[[i]] <- D[[i]]*Reduce("+", mapply(function(dil, Dl, Al) dil*(1-Dl)*Al, delta[[i]], D, A, SIMPLIFY=FALSE)[-i] )
    phi0x0 <- exp(V0)/x0
    G <- mapply(function(pi,Ei) phi0x0*pi - Ei, cost, E, SIMPLIFY=FALSE)
    
    # Calculate W (list) and its dimensions
    W <- mapply(function(vi, xi, gi, Gi) (vi - log(xi/gi + 1) - log(Gi)), V, X, gamma, G, SIMPLIFY=FALSE)
    nDrawsInter <- max(sapply(W, function(w) ifelse(is.array(w), dim(w)[2], 0) ))
    nDrawsIntra <- max(sapply(W, function(w) ifelse(is.array(w) && length(dim(w))==3, dim(w)[3], 0) ))
    
    # Calculate components of the jacobian for the whole sample 
    # J is a list of lists ( K x K ), where each element can be
    #   a vector, matrix or 3-dim array.
    J <- vector(mode="list", nAlt)
    for(i in 1:nAlt){
      J[[i]] <- vector(mode="list", nAlt)
      for(j in 1:nAlt){
        if(i==j) J[[i]][[j]] <- C[[i]] + (phi0x0*cost[[i]]^2/x0 + E[[i]])/G[[i]]
        if(i!=j) J[[i]][[j]] <- (phi0x0/x0*cost[[i]]*cost[[j]] - delta[[i]][[j]]*D[[i]]*D[[j]])/G[[i]]
      }
    }
    
    # Calculate determinant of each observation jacobian
    if(nDrawsInter==0 & nDrawsIntra==0){
      Jdet <- rep(1, nObs)
      for(n in which(Mv>0)){
        Jn <- matrix(0, Mv[n], Mv[n])
        ii <- which(Mm[n,])
        for(i in 1:Mv[n]) for(j in 1:Mv[n]) Jn[i,j] <- J[[ii[i]]][[ii[j]]][n]
        Jdet[n] <- det(Jn)
      }
    } else {
      Jdet <- array(1, dim=c(nObs, nDrawsInter+(nDrawsInter==0), 
                             nDrawsIntra+(nDrawsIntra==0)))
      for(n in which(Mv>0)){
        Jn <- matrix(0, Mv[n], Mv[n])
        ii <- which(Mm[n,])
        for(k in 1:dim(Jdet)[2]) for(l in 1:dim(Jdet)[3]){
          for(i in 1:Mv[n]) for(j in 1:Mv[n]){
            Jij <- J[[ii[i]]][[ii[j]]]
            if(is.array(Jij) && length(dim(Jij))==3) Jn[i,j] <- Jij[n, k, l]
            if(is.matrix(Jij)) Jn[i,j] <- Jij[n, k]
            if(is.vector(Jij)) Jn[i,j] <- Jij[n]
          }
          Jdet[n, k, l] <- det(Jn)
        }
      }
      if(nDrawsInter<=1 & nDrawsIntra<=1) Jdet <- as.vector(Jdet)
      if(nDrawsInter>1  & nDrawsIntra<=1) Jdet <- matrix(Jdet, 
                                                         nrow=nObs, 
                                                         ncol=nDrawsInter)
    }
    
    # Likelihood
    ll <- log(Jdet)
    ll <- ll + Reduce("+", mapply(function(w, m, a) dnorm(-w, sd=sigma, log=TRUE)*m*a, W, M, A, SIMPLIFY=FALSE))
    ll <- ll + Reduce("+", mapply(function(w, m, a) pnorm(-w, sd=sigma, log.p=TRUE)*(1-m)*a, W, M, A, SIMPLIFY=FALSE))
    
    return(exp(ll))
  }
  
  # ------------ #
  #### OUTPUT ####
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- apollo_emdc1(emdc_settings, functionality="estimate")
    
    ### Create report
    r <- c("Times available", "Times chosen", "% chosen when avail.", 
           "Avg. consump. when avail.", "Avg. consump. when chosen")
    M <- matrix(0, nrow=5, ncol=length(V), dimnames=list(r, names(X)))
    for(a in names(X)){
      # Expand avail
      if(length(A[[a]])==1) A[[a]] <- rep(A[[a]], length(X[[a]]))
      M[1, a] <- sum(A[[a]])            # Times available
      M[2, a] <- sum(X[[a]]>0)          # Times chosen
      M[3, a] <- M[2, a]/M[1, a]*100    # % chosen when avail
      M[4, a] <- mean(X[[a]][A[[a]]>0]) # Avg. consump. when avail.
      M[5, a] <- mean(X[[a]][X[[a]]>0]) # Avg. consump. when chosen.
    }
    content <- list(round(M,2))
    if(any(M[5,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
    if(any(M[3,]==100)) content[[length(content)+1]] <- "Warning: some alternatives are always chosen when available!"
    if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                               "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
    
    if(functionality=='report') ans <- list(like=ans, param=capture.output(print(M)))
    
    return(ans)
  }
  
  # ---------------- #
  #### PREDICTION ####
  # ---------------- #
  if(functionality=="prediction"){
    ### Find out dimensionality
    f <- function(x, maxD=c(-1,-1,-1)){
      d <- maxD
      if(is.vector(x)) d <- c(length(x), 1, 1)
      if(is.matrix(x)) d <- c(dim(x), 1)
      if(is.array(x) && length(dim(x))==3) d <- dim(x)
      maxD <- pmax(maxD, d)
      return(maxD)
    }
    maxD <- f(V0); maxD <- f(B, maxD); maxD <- f(A, maxD)
    for(k in 1:nAlt){
      maxD <- f(X[[k]], maxD)
      maxD <- f(cost[[k]], maxD)
      maxD <- f(V[[k]], maxD)
      maxD <- f(gamma[[k]], maxD)
      for(kk in 1:nAlt) maxD <- f(delta[[k]][[kk]], maxD)
    }; rm(f, k, kk)
    nInter <- maxD[2]
    nIntra <- maxD[3]
    rm(maxD)
    
    ### Function to extract required obs
    extractN <- function(y, n, inter, intra){
      if(is.list(y)) return(sapply(y, extractN, n=n, inter=inter, intra=intra))
      if(is.vector(y)){
        if(length(y)==1) return(y) else return(y[n])
      }
      if(is.matrix(y)) return(y[n,inter])
      if(is.array(y) && length(dim(y))==3) return(y[n,inter,intra])
      stop('Element [', n, ', ', inter, ', ', intra, '] does not exist.')
    }; environment(extractN) <- new.env(parent=environment())
    
    ### Make cluster
    # Each worker must contain phi0, V, gamma, A, delta, cost, B, nObs, nAlt, nInter, nIntra, timeLimit, extractN, & algoX
    test <- !is.null(apollo_inputs$apollo_control) && !is.null(apollo_inputs$apollo_control$nCores)
    if(test) nCores <- apollo_inputs$apollo_control$nCores else nCores <- 1
    # Split and copy into each cluster
    assignedCore <- rep(1:nCores, each=floor(nObs/nCores))
    if(length(assignedCore)<nObs) assignedCore <- c(assignedCore, rep(nCores, nObs-length(assignedCore)))
    inputPieceFile <- rep("", nCores)
    # Write pieces to disk
    for(i in 1:nCores){
      inputPieceFile[i] <- tempfile()
      r <- assignedCore==i
      L <- list(
        phi0 = exp(apollo_keepRows(V0, r=r)),
        V    = lapply(V    , apollo_keepRows, r=r), 
        gamma= lapply(gamma, apollo_keepRows, r=r), 
        A    = lapply(A    , apollo_keepRows, r=r), 
        delta= lapply(delta, function(d) lapply(d, apollo_keepRows, r=r)), 
        cost = lapply(cost , apollo_keepRows, r=r), 
        B    = apollo_keepRows(B, r=r), 
        nObs = nObs, 
        nAlt = nAlt, 
        nInter= nInter, 
        nIntra= nIntra, 
        timeLimit= timeLimit,
        extractN = extractN
      )
      wroteOK <- tryCatch({
        saveRDS(L, file=inputPieceFile[i])
        TRUE
      }, warning=function(w) FALSE, error=function(e) FALSE)
      if(!wroteOK) stop("Workers could not be set up for prediction")
    }
    # Create cluster
    cl <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(cl), add=TRUE)
    parallel::clusterEvalQ(cl, library('Rsolnp'))
    # Read pieces in each cluster and delete temporary files
    parallel::parLapply(cl, as.list(inputPieceFile), fun=function(fileName){
      tmp <- globalenv()
      if(!file.exists(fileName)) stop("A piece of data is missing from disk.")
      L <- tryCatch(readRDS(fileName), warning = function(w) FALSE, error = function(e) FALSE)
      if(is.logical(L) && !L) stop("A piece of data could not be loaded from disk.")
      for(i in 1:length(L)) assign(names(L)[i], L[[i]], envir=tmp)
    })
    rm(L)
    unlink(inputPieceFile)
    
    ### Generate random disturbances and split them across threads
    # Generate uniform draws,transform them to normal, and put them in an nObs x nAlt x nRep cube
    set.seed(24)
    #tmp <- apollo:::apollo_mlhs(nRep, nAlt, nObs)
    tmp <- apollo_mlhs(nRep, nAlt, nObs)
    U   <- array(0, dim=c(nObs, nAlt, nRep))
    for(iRep in 1:nRep) for(k in 1:nAlt) U[,k,iRep] <- qnorm(tmp[((iRep-1)*nObs+1):(iRep*nObs), k], sd=sigma)
    rm(tmp)
    # Split the cube of draws into a list with three levels: [[rep]][[core]][[alt]]
    EPS <- vector(mode='list', length=nRep)
    for(iRep in 1:nRep){
      EPS[[iRep]] <- list()
      for(i in 1:nCores){
        r <- assignedCore==i
        EPS[[iRep]][[i]] <- list()
        for(k in 1:nAlt) EPS[[iRep]][[i]][[k]] <- U[r,k,iRep]
      }
    }; rm(U)
    
    ### Repetition loop
    # Initialise output variables
    alts  <- c('outside', names(X))
    Xmean <- matrix(0, nrow=nObs, ncol=nAlt+1)
    Mmean <- Xmean; Xsd <- Xmean; Msd <- Xmean
    if(!apollo_inputs$silent){ iHalf <- ceiling(nRep/2); iDot <- ceiling(nRep/10) }
    for(iRep in 1:nRep){
      if(!apollo_inputs$silent && iRep==1) cat(" 0%")
      # Forecast
      parallel::parLapply(cl, EPS[[iRep]], fun=function(eps){tmp <- globalenv(); assign('eps', eps, envir=tmp)})
      X1 <- parallel::clusterEvalQ(cl, {
        # Load variables from global environment to prevent NOTES
        eps  <- get("eps" , envir=globalenv(), inherits=FALSE)
        phi0 <- get("phi0", envir=globalenv(), inherits=FALSE)
        # Update nObs
        nObs <- max(sapply(eps, length))
        # Calculate base utility
        phi <- mapply(function(Vi, ei) exp(Vi + ei), V, eps, SIMPLIFY=FALSE)
        # Initialise output
        xO <- array(NA, dim=c(nObs, nInter, nIntra))
        X  <- list()
        for(k in 1:nAlt) X[[k]] <- xO
        ### Loop over draws and observations
        for(inter in 1:nInter) for(intra in 1:nIntra) for(n in 1:nObs) {
          # Extract values
          ph0 <- extractN(phi0 , n, inter, intra)
          ph  <- extractN(phi  , n, inter, intra)
          g   <- extractN(gamma, n, inter, intra)
          p   <- extractN(cost , n, inter, intra)
          b   <- extractN(B    , n, inter, intra)
          d   <- matrix(0, nAlt, nAlt)
          for(j in 1:(nAlt-1)) for(i in (j+1):nAlt) d[i,j] <- extractN(delta[[i]][[j]], n, inter, intra)
          a   <- extractN(A, n, inter, intra)>0 # make sure it is logical
          # Remove unavailable alternatives
          if(!all(a)){
            ph <- ph[a]
            g  <- g[a]
            p  <- p[a]
            d  <- d[a,a]
          }
          # Starting solution proportional to phi/p
          x <- c(ph0, ph/p)
          x <- x/sum(x)*b/c(1,p)
          # Set time limit
          setTimeLimit(cpu=timeLimit, elapsed=timeLimit, transient=TRUE)
          # Optimisation with time limit
          eMsg <- "reached elapsed time limit|reached CPU time limit"
          tryCatch({
            m <- Rsolnp::solnp(pars=x, 
                               fun=function(x){ D <- 1 - exp(-x[-1])
                               ans <- ph0*log(x[1]) + sum( ph*g*log(x[-1]/g + 1) ) + sum((d*D)%*%D)
                               return(-ans)}, 
                               eqfun=function(x) x[1] + sum(x[-1]*p), 
                               eqB=b, LB=rep(0, length(x)), 
                               control=list(trace=0, outer.iter=200, inner.iter=400))
            x <- m$par
          }, error = function(e) if(grepl(eMsg, e$message)) x <- rep(NA, length(x)) else stop(e) )
          # Remove time limit
          setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE)
          # Expand x if there were some unavailable alternatives
          if(!all(a)){
            tmp <- rep(0, length(a)+1)
            tmp[c(TRUE, a)] <- x
            x   <- tmp
            rm(tmp)
          }
          # Store results
          xO[n, inter, intra] <- x[1]
          for(k in 1:nAlt) X[[k]][n, inter, intra] <- x[1+k]
        }
        ### Return
        xO <- apply(xO, MARGIN=1, mean)
        X <- sapply(X, function(xk) apply(xk, MARGIN=1, mean))
        X <- cbind(xO, X)
        colnames(X) <- c('outside', names(phi))
        X
      })
      X1 <- do.call(rbind, X1)
      # Store results
      Xmean <- Xmean + X1/nRep
      Mmean <- Mmean + (X1>0)/nRep
      Xsd   <- Xsd   + apply(X1, MARGIN=2, function(x) ((x-mean(x))^2)/nRep)
      Msd   <- Msd   + apply(X1, MARGIN=2, function(x){ x <- x>0; ((x-mean(x))^2)/nRep} )
      # Print progress
      if(!apollo_inputs$silent){
        if(iRep==nRep) cat("100%") else if(iRep%%iHalf==0) cat("50%") else if(iRep%%iDot==0) cat(".")
      }
    }
    
    ### Prepare output
    Xsd <- apply(Xsd, MARGIN=2, sqrt)
    Msd <- apply(Msd, MARGIN=2, sqrt)
    return( cbind(Xmean, Xsd, Mmean, Msd) )
  }
  
  # End of function
  stop("Invalid value of argument 'functionality'")
}