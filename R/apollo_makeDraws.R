#' Create draws for models with mixing
#'
#' Creates a list containing all draws necessary to estimate a model with mixing.
#'
#' This function creates a list whose elements are the sets of draws requested by the user for use in a model with mixing.
#' If the model does not include mixing, then it is not necessary to run this function.
#' The number of draws have a massive impact on memory usage and estimation time. Memory usage and number of computations
#' scale geometrically as N*inter_nDraws*intra_nDraws (where N is the number of observations). Special care should be taken
#' when using both inter and intra draws, as memory usage can easily reach the GB order of magnitude. Also, keep in
#' mind that using several threads (i.e. multicore) at least doubles the memory usage.
#' This function returns a list, with each element representing a random component of the mixing model (except the last
#' element, which is a copy of the argument \code{apollo_draws}). The dimensions of the array depends on the type of draws used.
#' \enumerate{
#'            \item If only inter-individual draws are used, then draws are stored as 2-dimensional arrays (i.e. matrices).
#'            \item If intra-individual draws are used, then draws are stored as 3-dimensional arrays.
#'            \item The first dimension of the arrays (rows) correspond with the observations in the database.
#'            \item The second dimension of the arrays (columns) correspond to the number of inter-individual draws.
#'            \item The third dimension of the arrays correspond to the number of intra-individual draws.
#' }
#' @param apollo_control List. Contains options for the estimation
#'                    See \link{apollo_validatecontrol} for details.
#' @param apollo_draws List of arguments describing the inter and intra individual draws.
#'                  \describe{
#'                    \item{inter_drawsType}{Character. Type of inter-individual draws ('MLHS', 'halton' or 'pmc').}
#'                    \item{inter_nDraws}{Numeric scalar (>=0). Number of inter-individual draws per individual. Set to 0 if not using them.}
#'                    \item{inter_unifDraws}{Character vector. Names of uniform-distributed inter-individual draws.}
#'                    \item{inter_normDraws}{Character vector. Names of normaly distributed inter-individual draws.}
#'                    \item{intra_drawsType}{Character. Type of intra-individual draws ('MLHS', 'halton' or 'pmc').}
#'                    \item{intra_nDraws}{Numeric scalar (>=0). Number of intra-individual draws per individual. Set to 0 if not using them.}
#'                    \item{intra_unifDraws}{Character vector. Names of uniform-distributed intra-individual draws.}
#'                    \item{intra_normDraws}{Character vector. Names of normaly distributed intra-individual draws.}
#'                  }
#' @param database data.frame. Model data.
#' @param silent Boolean. If true, then no information is printed to console or default output. FALSE by default.
#' @return List. Each element is an array of draws representing a random component of the mixing model. The last element is a copy of argument \code{apollo_draws}.
apollo_makeDraws=function(apollo_control, apollo_draws, database, silent=FALSE){
  if(!apollo_control$mixing){
    warning("No need to call apollo_makeDraws if no mixing is used.")
    return(NA)
  }

  # ################################## #
  #### Validation of apollo_draws      ####
  # ################################## #
  interComplete <- TRUE
  if( is.null(apollo_draws$inter_drawsType) || !(apollo_draws$inter_drawsType %in% c('halton','mlhs','pmc')) ) interComplete <- FALSE
  if( is.null(apollo_draws$inter_nDraws) || !is.numeric(apollo_draws$inter_nDraws) || apollo_draws$inter_nDraws<0 ) interComplete <- FALSE
  if( (is.null(apollo_draws$inter_unifDraws) || !is.character(apollo_draws$inter_unifDraws) || length(apollo_draws$inter_unifDraws)==0) &
      (is.null(apollo_draws$inter_normDraws) || !is.character(apollo_draws$inter_normDraws) || length(apollo_draws$inter_normDraws)==0) ) interComplete <- FALSE
  intraComplete <- TRUE
  if( is.null(apollo_draws$intra_drawsType) || !(apollo_draws$intra_drawsType %in% c('halton','mlhs','pmc')) ) intraComplete <- FALSE
  if( is.null(apollo_draws$intra_nDraws) || !is.numeric(apollo_draws$intra_nDraws) || apollo_draws$intra_nDraws<0 ) intraComplete <- FALSE
  if( (is.null(apollo_draws$intra_unifDraws) || !is.character(apollo_draws$intra_unifDraws) || length(apollo_draws$intra_unifDraws)==0) &
      (is.null(apollo_draws$intra_normDraws) || !is.character(apollo_draws$intra_normDraws) || length(apollo_draws$intra_normDraws)==0) ) intraComplete <- FALSE
  if(!interComplete & !intraComplete){
    stop("Invalid draws settings in apollo_draws. See ?apollo_makeDraws.")
  }
  if(!interComplete){
    apollo_draws$inter_drawsType="mlhs"
    apollo_draws$inter_nDraws=0
    apollo_draws$inter_unifDraws=c()
    apollo_draws$inter_normDraws=c()
    if(!silent) cat("No inter-individual draws created.\n")
  }
  if(!intraComplete){
    apollo_draws$intra_drawsType="mlhs"
    apollo_draws$intra_nDraws=0
    apollo_draws$intra_unifDraws=c()
    apollo_draws$intra_normDraws=c()
    if(!silent) cat("No intra-individual draws created.\n")
  }



  # ################################## #
  #### Initialistion                ####
  # ################################## #

  panelData <- apollo_control$panelData
  indivID   <- database[,apollo_control$indivID]

  if(is.null(apollo_control$seed_draws)) apollo_control$seed_draws=13
  set.seed(apollo_control$seed_draws)

  d <- apollo_draws
  nObs <- length(indivID)
  if(!panelData) indivID <- 1:nObs
  nIndiv <- length(unique(indivID))

  namesInter <- c(d$inter_unifDraws, d$inter_normDraws)
  namesIntra <- c(d$intra_unifDraws, d$intra_normDraws)
  dimInter <- length(namesInter)
  dimIntra <- length(namesIntra)

  d$inter_drawsType <- tolower(d$inter_drawsType)
  d$intra_drawsType <- tolower(d$intra_drawsType)
  invalidInter <- ( dimInter>0 & !(d$inter_drawsType %in% c('halton','mlhs','pmc')) )
  invalidIntra <- ( dimIntra>0 & !(d$intra_drawsType %in% c('halton','mlhs','pmc')) )
  if( invalidInter | invalidIntra) stop('Invalid type of draws. Use "halton", "mlhs" or "pmc".')

  if(d$inter_nDraws==0 | dimInter==0) {d$inter_nDraws <- 1; dimInter=0}
  if(d$intra_nDraws==0 | dimIntra==0) {d$intra_nDraws <- 1; dimIntra=0}

  drawsList <- list()
  if(!silent) cat('Creating draws ')

  # ################################## #
  #### Inter-individual draws       ####
  # ################################## #

  if(dimInter>0){
    if(apollo_draws$inter_drawsType=='halton') draws <- randtoolbox::halton(d$inter_nDraws*nIndiv,
                                                                         dimInter)
    if(apollo_draws$inter_drawsType=='mlhs') draws <- apollo_mlhs(d$inter_nDraws,dimInter,nIndiv)
    if(apollo_draws$inter_drawsType=='pmc') draws <- matrix(stats::runif(nIndiv*d$inter_nDraws*dimInter),
                                                         nrow=nIndiv*d$inter_nDraws, ncol=dimInter,
                                                         byrow=TRUE)
    draws <- as.matrix(draws)
    colnames(draws) <- c(d$inter_unifDraws, d$inter_normDraws)

    if(length(d$inter_normDraws)>0){
      for(i in (length(d$inter_unifDraws)+1):dimInter) draws[,i] <- stats::qnorm(draws[,i])
    }

    obsPerIndiv <- as.vector(table(indivID))
    for(d1 in 1:dimInter){

      M <- matrix(0, nrow=nObs, ncol=d$inter_nDraws)
      row1 <- 1
      for(i in 1:nIndiv){
        row2 <- row1 + obsPerIndiv[i] - 1
        M[row1:row2,] <- matrix(draws[((i-1)*d$inter_nDraws+1):(i*d$inter_nDraws),d1],
                                nrow=row2-row1+1, ncol=d$inter_nDraws, byrow=TRUE)
        row1 <- row2 + 1
      }
      for(d2 in 1:dimInter){
        C <- array(0, dim=c(nObs, d$inter_nDraws, d$intra_nDraws))
        for(j in 1:d$intra_nDraws) C[,,j] <- M
      }
      drawsList[[namesInter[d1]]] <- C
      if(!silent) cat('.')
    }
  }

  # ################################## #
  #### Intra-individual draws       ####
  # ################################## #

  if(dimIntra>0){

    if(d$intra_drawsType=='halton'){
      if(dimInter>0 & apollo_draws$inter_drawsType=='halton'){
        draws <- randtoolbox::halton(d$intra_nDraws*nObs,dimInter+dimIntra)
        draws <- draws[,(dimInter+1):(dimInter+dimIntra)]
      } else {
        draws <- randtoolbox::halton(d$intra_nDraws*nObs,dimIntra) # If inter is not halton
      }
      draws <- as.matrix(draws)
    }
    if(d$intra_drawsType=='mlhs'){
      draws <- apollo_mlhs(d$intra_nDraws,dimIntra,nObs)
    }
    if(d$intra_drawsType=='pmc'){
      draws <- stats::runif(d$intra_nDraws*nObs*dimIntra)
      draws <- matrix(draws, nrow=d$intra_nDraws*nObs, ncol=dimIntra)
    }

    if(length(d$intra_normDraws)>0){
      for(d2 in (length(d$intra_unifDraws)+1):dimIntra) draws[,d2] <- stats::qnorm(draws[,d2])
    }

    for(d2 in 1:dimIntra){
      C <- array(0, dim=c(nObs, d$inter_nDraws, d$intra_nDraws))
      for(n in 1:nObs){
        depthInterval <- ((n-1)*d$intra_nDraws+1):(n*d$intra_nDraws)
        C[n,1,] <- draws[depthInterval,d2]
      }
      if(d$inter_nDraws>1) for(dinter in 2:d$inter_nDraws) C[,dinter,] <- C[,1,]
      drawsList[[namesIntra[d2]]] <- C
      if(!silent) cat('.')
    }

  }


  # ################################## #
  #### Returning draws              ####
  # ################################## #

  if(!silent) cat(' Done\n')
  if(!panelData & dimInter>0){
    warning('Inter-person draws are being used without a panel structure.')
  }

  if(d$intra_nDraws<2){
    j=1
    while(j<=length(drawsList)){
      drawsList[[j]] <- colSums(aperm(drawsList[[j]], perm=c(3,1,2)))/dim(drawsList[[j]])[3]
      j=j+1
    }
  }

  drawsList[['apollo_draws']] <- apollo_draws

  if(!silent & !interComplete & intraComplete) cat("Intra-individual draws created.")
  if(!silent & interComplete & !intraComplete) cat("Inter-individual draws created.")
  if(!silent & interComplete & intraComplete) cat("Inter and intra-individual draws created.")


  return(drawsList)
}
