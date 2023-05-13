#' Extended MDC
#' 
#' Calculates the likelihood function of the extended MDC model. Can also predict and validate inputs.
#' 
#' This model extends the traditional multiple discrete-continuous (MDC) framework by (i) dropping
#' the need to define a budget, (ii) making the marginal utility of the outside good deterministic, 
#' and (iii) including complementarity and substitution in the model formulation. See the following 
#' working paper for more details:
#' 
#' Palma, D. & Hess, S. (Working Paper) Some adaptations of Multiple Discrete-Continuous 
#' Extreme Value (MDCEV) models for a computationally tractable treatment of complementarity 
#' and substitution effects, and reduced influence of budget assumptions
#' 
#' Avilable at: http://stephanehess.me.uk/publications.html
#' 
#' @param emdc_settings List of settings for the model. It includes the following.
#'                        \itemize{
#'                          \item \strong{\code{continuousChoice}}: Named list of numeric vectors. Amount consumed of each inside good. Outside good must not be included. Can also be called "X".
#'                          \item \strong{\code{avail}}: Named list of numeric vectors. Availability of each product. Can also be called "A".
#'                          \item \strong{\code{utilityOutside}}: Numeric vector (or matrix or array). Shadow price of the budget. Must be normalised to 0 for at least one individual. Default is 0 for every observation. Can also be called "V0".
#'                          \item \strong{\code{utilities}}: Named list of numeric vectors (or matrices or arrays). Base utility of each product. Can also be called "V".
#'                          \item \strong{\code{gamma}}: Named list of numeric vectors. Satiation parameter of each product.
#'                          \item \strong{\code{sigma}}: Numeric scalar. Scale parameter.
#'                          \item \strong{\code{delta}}: Lower triangular numeric matrix, or list of lists. Complementarity/substitution parameter.
#'                          \item \strong{\code{cost}}: Named list of numeric vectors. Price of each product.
#'                          \item \strong{\code{nRep}}: Scalar positive integer. Number of repetitions used when predictiong
#'                          \item \strong{\code{nIter}}: Vector of two positive integers. Number of maximum iterations used during prediction, for the upper and lower iterative levels.
#'                          \item \strong{\code{tolerance}}: Positive scalar Tolerance of the prediction algorithm.
#'                          \item \strong{\code{rawPrediction}}: Scalar logical. When functionality is equal to "prediction", it returns the full set of simulations. Defaults is FALSE.
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
apollo_emdc2 <- function(emdc_settings, functionality="estimate"){
  # Rename input if necessary
  map <- c(X = "continuousChoice", A = "avail", 
           V0= "utilityOutside",   V = "utilities")
  for(i in 1:length(map)) if(!is.null(emdc_settings[[map[i]]])){
    emdc_settings[[names(map)[i]]] <- emdc_settings[[map[i]]]
    emdc_settings[[map[i]]]        <- NULL
  }; rm(i, map)
  
  # Check input
  mandatory <- c("X", "A", "V0", "V", "gamma", "delta", "cost")
  optional  <- list(sigma=1, nRep=50, nIter=c(50,50), tolerance=0.1, rawPrediction=FALSE)
  if(!all(mandatory %in% names(emdc_settings))) stop("Some mandatory elements are missing from 'emdc_settings'.")
  for(i in names(optional)) if(!(i %in% names(emdc_settings))) emdc_settings[[i]] <- optional[[i]]
  test <- is.list(emdc_settings$X) && all(sapply(emdc_settings$X, is.vector))
  if(!test) stop("'X' should be a list of numeric vectors.")
  test <- is.list(emdc_settings$A) && all(sapply(emdc_settings$A, is.vector))
  if(!test) stop("'A' should be a list of numeric vectors.")
  test <- is.vector(emdc_settings$V0) | is.array(emdc_settings$V0) 
  if(!test) stop("'V0' should be a numeric vector, matrix or array.")
  test <- is.list(emdc_settings$V) && all(sapply(emdc_settings$V, function(v) is.vector(v) | is.array(v)))
  if(!test) stop("'V' should be a list of numeric vectors, matrices or arrays.")
  test <- is.list(emdc_settings$gamma) && all(sapply(emdc_settings$gamma, function(g) is.vector(g) | is.array(g)))
  if(!test) stop("'gamma' should be a list of numeric vectors, matrices or arrays.")
  test <- is.vector(emdc_settings$sigma) && length(emdc_settings$sigma)==1
  if(!test) stop("'sigma' should be a numeric scalar.")
  test <- is.list(emdc_settings$cost) && all(sapply(emdc_settings$cost, is.vector))
  if(!test) stop("'cost' should be a list of numeric vectors.")
  test <- is.vector(emdc_settings$nRep) && length(emdc_settings$nRep)==1 && emdc_settings$nRep>0
  if(!test) stop("'nRep' should be a positive scalar integer.")
  test <- is.vector(emdc_settings$nIter) && length(emdc_settings$nIter)==2 && all(emdc_settings$nIter>0)
  if(!test) stop("'nIter' should be a numeric vector of two positive elements.")
  test <- is.vector(emdc_settings$tolerance) && length(emdc_settings$tolerance)==1 && emdc_settings$tolerance>0
  if(!test) stop("'tolerance' should be a positive numeric scalar.")
  test <- is.vector(emdc_settings$rawPrediction) && length(emdc_settings$rawPrediction)==1 && is.logical(emdc_settings$rawPrediction)
  if(!test) stop("'rawPrediction should be a scalar logical (TRUE OR FALSE).")
  rm(mandatory, optional, test, i)
  
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE, debug=FALSE))
  test <- !is.null(apollo_inputs$silent)
  if(test) silent <- apollo_inputs$silent else silent <- FALSE
  test <- !is.null(apollo_inputs$debug )
  if(test) debug <- apollo_inputs$debug else debug <- FALSE
  
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
  nRep  <- nRep  # this prevents a NOTE when checking the package.
  nIter <- nIter # this prevents a NOTE when checking the package.
  sigma <- sigma # this prevents a NOTE when checking the package.
  tolerance <- tolerance # this prevents a NOTE when checking the package.
  rawPrediction <- rawPrediction # this prevents a NOTE when checking the package.
  A     <- A[names(X)]
  V     <- V[names(X)]
  gamma <- gamma[names(X)]
  cost  <- cost[names(X)]
  nAlt  <- length(X)
  nObs  <- c(sapply(X, length), sapply(A, length), 
             sapply(gamma, length), sapply(cost, length))
  nObs  <- max(nObs)
  
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
    rm(d)
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
    if(!is.vector(nIter) || length(nIter)!=2 || any(nIter<=0)) stop("nIter must be a vector of length 2 with positive integers")
    if(!is.vector(tolerance) || length(tolerance)!=1 || tolerance<0) stop("tolerance must be a non-negative scalar")
    if(min(sapply(X, length))==0) stop("At least one element in X has length 0")
    
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
    #apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    #apollo:::apollo_addLog("Overview of choices for eMDCEV model component:", content, apolloLog)
    
    testL <- apollo_emdc2(emdc_settings, functionality="estimate")
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
    # Create useful variables
    M  <- lapply(X, function(x) x>0)
    Mm <- simplify2array(M)
    Mv <- rowSums(Mm)
    
    # Calculate C, D, E, G
    C <- mapply(function(xi, gi) 1/(xi + gi), X, gamma, SIMPLIFY=FALSE)
    D <- lapply(X, function(xi) exp(-xi))
    E <- vector(mode="list", nAlt)
    for(i in 1:nAlt) E[[i]] <- D[[i]]*Reduce("+", mapply(function(dil, Dl, Ai) dil*(1-Dl)*Ai, delta[[i]], D, A, SIMPLIFY=FALSE)[-i] )
    eV0 <- exp(V0)
    G <- mapply(function(pi,Ei) eV0*pi - Ei, cost, E, SIMPLIFY=FALSE)
    
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
        if(i==j) J[[i]][[j]] <- C[[i]] + E[[i]]/G[[i]]
        if(i!=j) J[[i]][[j]] <- -delta[[i]][[j]]*D[[i]]*D[[j]]/G[[i]]
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
    ans <- apollo_emdc2(emdc_settings, functionality="estimate")
    
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
    #apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
    #apollo:::apollo_addLog("Overview of choices for eMDCEV model component:", content, apolloLog)
    
    return(ans)
  }
  
  # ---------------- #
  #### PREDICTION ####
  # ---------------- #
  if(functionality=="prediction"){
    #cat("\nNow producing forecasts from eMDCEV model component.")
    #cat("\nFor each product (columns) and each observation (rows),\n the following will be returned in a matrix:")
    #cat("\n 1. Expected consumption")
    #cat("\n 2. Std. dev. of expected consumption")
    #cat("\n 3. Expected probability of discrete consumptions")
    #cat("\n 4. Std. dev. of exp. prob. of disc. consumptions\n")
    
    ### Fixed point forecasting for X
    fixedPointX <- function(phi0, phi, gamma, delta, Xstart, cost, A, maxIter, tol){
      nAlt2<- length(gamma)
      i1   <- 1
      rmse <- Inf
      X0   <- Xstart
      X1   <- Xstart
      if(nAlt2==0) return(list()) # if all alternatives are fixed to zero consumption
      while(i1<=maxIter[1] && rmse>tol){ # if there is at least one alternative
        for(j in 1:nAlt2){
          i2 <- 1
          D  <- mapply(function(xl, djl, al) djl*(1-exp(-xl))*al, X1, delta[[j]], A, SIMPLIFY=FALSE)
          Dj <- Reduce("+", D[-j])
          tmp<- X0[[j]]
          err<- Inf
          while(i2<=maxIter[2] && sqrt(mean(err))>tol){
            Ej      <- exp(-X1[[j]])*Dj
            X1[[j]] <- gamma[[j]]*(phi[[j]]/(phi0*cost[[j]] - Ej) - 1)*A[[j]]
            err     <- (X1[[j]] - tmp)^2
            if(is.array(err)){
              if(is.matrix(err)) err <- rowMeans(err)
              if(length(dim(err))==3) err <- apply(err, MARGIN=c(1,2), mean)
            }
            tmp     <- X1[[j]]
            i2      <- i2 + 1
          }
          duErr <- phi[[j]]/(X1[[j]]/gamma[[j]] + 1) + Ej - phi0*cost[[j]] # dUj - phi0*costj should be 0 for consumed ones
          X1[[j]][ X1[[j]]<0 | abs(duErr)>tol | sqrt(err)>tol ] <- 0
        }
        rmse <- Reduce("+", mapply(function(x1, x0) (x1-x0)^2, X1, X0, SIMPLIFY=FALSE))
        if(is.array(rmse)){
          if(length(dim(rmse))==3) rmse <- apply(rmse, MARGIN=c(1,2), mean)
          if(is.matrix(rmse)) rmse <- mean(sqrt(colMeans(rmse)))
        } else rmse <- sqrt(mean(rmse))
        X0   <- X1
        i1   <- i1 + 1
      }
      rmseSum <<- rmseSum + rmse
      names(X0) <- names(Xstart)
      X0 <- lapply(X0, function(xk) if(is.array(xk)) apply(xk, MARGIN=1, mean) else xk)
      return(X0)
    }
    
    ### Initialise output variables
    alts  <- names(X)
    if(rawPrediction) Xmean <- array(NA, dim=c(nObs, nAlt, nRep), dimnames=list(NULL, alts, NULL)) else {
      Xmean <- list()
      Mmean <- list()
      Xsd   <- list()
      Msd   <- list()
      for(i in 1:nAlt){
        Xmean[[paste0(alts[i], "_ECont")]]  <- rep(0, nObs)
        Mmean[[paste0(alts[i], "_EDisc")]]  <- rep(0, nObs)
        Xsd  [[paste0(alts[i], "_sdCont")]] <- rep(0, nObs)
        Msd  [[paste0(alts[i], "_sdDisc")]] <- rep(0, nObs)
      }
    }
    rmseSum <- 0
    
    ### Repetition loop
    set.seed(24)
    phi0 <- exp(V0)
    # Generate random disturbances
    # eps is a list with nRep elements, each of which is a list with nAlt elements
    #tmp <- apollo:::apollo_mlhs(nRep, nAlt, nObs)
    tmp <- apollo_mlhs(nRep, nAlt, nObs)
    tmp <- qnorm(tmp, sd=sigma)
    eps <- list()
    for(i in 1:nRep){
      eps[[i]] <- tmp[ (0:(nObs-1))*nRep + i, ]
      eps[[i]] <- split(eps[[i]], rep(1:nAlt, each=nObs))
    }; rm(tmp)
    for(r in 1:nRep){
      # Generate starting values
      if(r==1) X <- lapply(X, function(x) x*0) else {
        # Exponential
        L <- lapply(X1, mean)
        X <- lapply(L, function(l) -log(runif(nObs))/ifelse(l==0, 1, l) )
      }
      # Calculate base utilities and forecast
      phi <- mapply(function(Vi, ei) exp(Vi + ei), V, eps[[r]], SIMPLIFY=FALSE)
      X1  <- fixedPointX(phi0, phi, gamma, delta, X, cost, A, maxIter=nIter+r%%2, tol=tolerance) # adds an extra iter to odd repetitions.
      # Store results
      if(rawPrediction) Xmean[,,r] <- do.call(cbind, X1) else {
        Xmean <- mapply(function(x0,  x) x0  +                       x/nRep, Xmean, X1, SIMPLIFY=FALSE)
        Mmean <- mapply(function(m0,  x) m0  +                   (x>0)/nRep, Mmean, X1, SIMPLIFY=FALSE)
        Xsd   <- mapply(function(xsd, x) xsd +         ((x-mean(x))^2)/nRep,   Xsd, X1, SIMPLIFY=FALSE)
        Msd   <- mapply(function(msd, x) msd + (((x>0)-mean((x>0)))^2)/nRep,   Msd, X1, SIMPLIFY=FALSE)
      }
      if(!silent){ if(r==1) cat(" 0%") else if(r==nRep) cat("100%") else if(r%%ceiling(nRep/2)==0) cat("50%") else if(r%%ceiling(nRep/10)==0) cat(".") }
    }
    if(debug) cat(" conv_RMSE=", round(rmseSum/nRep,3), sep="")
    if(!silent) cat("\n")
    
    ### Prepare output
    if(rawPrediction) out <- Xmean else {
      Xsd   <- lapply(Xsd, sqrt)
      Msd   <- lapply(Msd, sqrt)
      Xmean <- do.call(cbind, Xmean)
      Mmean <- do.call(cbind, Mmean)
      Xsd   <- do.call(cbind,   Xsd)
      Msd   <- do.call(cbind,   Msd)
      out   <- cbind(Xmean, Xsd, Mmean, Msd)
    }
    return(out)
  }
  
  # End of function
  stop("Invalid value of argument 'functionality'")
}