#' Creates draws for models with mixing
#'
#' Creates a list containing all draws necessary to estimate a model with mixing.
#'
#' Internal use only. Called by \code{apollo_validateInputs}.
#' #' This function creates a list whose elements are the sets of draws requested by the user for use in a model with mixing.
#' If the model does not include mixing, then it is not necessary to run this function.
#' The number of draws have a massive impact on memory usage and estimation time. Memory usage and number of computations
#' scale geometrically as N*interNDraws*intraNDraws (where N is the number of observations). Special care should be taken
#' when using both inter and intra draws, as memory usage can easily reach the GB order of magnitude. Also, keep in
#' mind that using several threads (i.e. multicore) at least doubles the memory usage.
#' This function returns a list, with each element representing a random component of the mixing model. The dimensions 
#' of the array depend on the type of draws used.
#' \enumerate{
#'            \item If only inter-individual draws are used, then draws are stored as 2-dimensional arrays (i.e. matrices).
#'            \item If intra-individual draws are used, then draws are stored as 3-dimensional arrays.
#'            \item The first dimension of the arrays (rows) correspond with the observations in the database.
#'            \item The second dimension of the arrays (columns) correspond to the number of inter-individual draws.
#'            \item The third dimension of the arrays correspond to the number of intra-individual draws.
#' }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If true, then no information is printed to console or default output. FALSE by default.
#' @return List. Each element is an array of draws representing a random component of the mixing model.
#' @importFrom randtoolbox halton sobol
apollo_makeDraws=function(apollo_inputs, silent=FALSE){
  apollo_control = apollo_inputs[["apollo_control"]]
  d              = apollo_inputs[["apollo_draws"]]
  database       = apollo_inputs[["database"]]

  if(!apollo_control$mixing) return(NA)

  # ################################## #
  #### Input validation             ####
  # ################################## #

  testEGen <- rep(FALSE, 5)
  testEUsr    <- !is.null(d$interDrawsType) && d$interDrawsType!="" && length(d$interDrawsType)==1 && ( !(tolower(d$interDrawsType) %in% c('halton','mlhs','pmc','sobol','sobolowen','sobolfauretezuka','sobolowenfauretezuka')) && exists(d$interDrawsType, envir=globalenv()) )
  testEGen[1] <- !is.null(d$interDrawsType) && d$interDrawsType!="" && length(d$interDrawsType)==1 && (tolower(d$interDrawsType) %in% c('halton','mlhs','pmc','sobol','sobolowen','sobolfauretezuka','sobolowenfauretezuka'))
  testEGen[2] <- !is.null(d$interNDraws) && is.numeric(d$interNDraws) && d$interNDraws>0
  testEGen[3] <- length(d$interUnifDraws)==0 || is.character(d$interUnifDraws)
  testEGen[4] <- length(d$interNormDraws)==0 || is.character(d$interNormDraws)
  testEGen[5] <- length(c(d$interUnifDraws, d$interNormDraws)) > 0

  testAGen <- rep(FALSE, 5)
  testAUsr    <- !is.null(d$intraDrawsType) && d$intraDrawsType!="" && length(d$intraDrawsType)==1 && ( !(tolower(d$intraDrawsType) %in% c('halton','mlhs','pmc','sobol','sobolowen','sobolfauretezuka','sobolowenfauretezuka')) && exists(d$intraDrawsType, envir=globalenv()) )
  testAGen[1] <- !is.null(d$intraDrawsType) && d$intraDrawsType!="" && length(d$intraDrawsType)==1 && (tolower(d$intraDrawsType) %in% c('halton','mlhs','pmc','sobol','sobolowen','sobolfauretezuka','sobolowenfauretezuka'))
  testAGen[2] <- !is.null(d$intraNDraws) && is.numeric(d$intraNDraws) && d$intraNDraws>0
  testAGen[3] <- length(d$intraUnifDraws)==0 || is.character(d$intraUnifDraws)
  testAGen[4] <- length(d$intraNormDraws)==0 || is.character(d$intraNormDraws)
  testAGen[5] <- length(c(d$intraUnifDraws, d$intraNormDraws)) > 0

  if(!testEUsr & !all(testEGen)){
    d$interDrawsType=NA
    d$interNDraws=0
    d$interUnifDraws=c()
    d$interNormDraws=c()
  }
  if(!testAUsr & !all(testAGen)){
    d$intraDrawsType=NA
    d$intraNDraws=0
    d$intraUnifDraws=c()
    d$intraNormDraws=c()
  }
  
  
  # Validate input
  if(!testEUsr & !all(testEGen) & !testAUsr & !all(testAGen)){
    if(!testEGen[1] | !testAGen[1]) stop("Type of draws must be 'halton', 'mlhs', 'pmc', 'sobol','sobolOwen','sobolFaureTezuka' or 'sobolOwenFaureTezuka' to generate draws, or the name of a variable containing user generated draws.")
    if(!testEGen[2] | !testAGen[2]) stop("Number of draws must be larger than 0.")
    if(!testEGen[5] | !testAGen[5]) stop("No names for the draws were specified.")
  }


  # ################################## #
  #### Initialisation                ####
  # ################################## #

  panelData <- apollo_control$panelData
  indivID   <- database[,apollo_control$indivID]

  if(is.null(apollo_control$seed_draws)) apollo_control$seed_draws=13
  set.seed(apollo_control$seed_draws)

  nObs <- length(indivID)
  if(!panelData) indivID <- 1:nObs
  nIndiv <- length(unique(indivID))
  obsPerIndiv <- as.vector(table(indivID))
  nInter <- d$interNDraws
  nIntra <- d$intraNDraws
  namesInter <- c(d$interUnifDraws, d$interNormDraws)
  namesIntra <- c(d$intraUnifDraws, d$intraNormDraws)
  dimInter <- length(namesInter)
  dimIntra <- length(namesIntra)
  if(!testEUsr) d$interDrawsType <- tolower(d$interDrawsType)
  if(!testAUsr) d$intraDrawsType <- tolower(d$intraDrawsType)

  # Function that expand the dimensions of a matrix of draws to whatever is necessary
  expand <- function(M, isInter){
    # Repeat rows if necesary
    if(isInter & nrow(M)<nObs){
      M1 <- matrix(0, nrow=nObs, ncol=ncol(M))
      r1 <- 1
      for(i in 1:nIndiv){
        r2 <- r1 + obsPerIndiv[i] - 1
        M1[r1:r2,] <- matrix(as.vector(M[i,]), nrow=r2-r1+1, ncol=ncol(M), byrow=TRUE)
        r1 <- r2 + 1
      }
      M <- M1
      rm(M1, r1, r2, i)
    }

    # Turn into cube if necessary
    if(nIntra>0){
      C <- array(0, dim=c(nObs, max(1,nInter), nIntra))
      if(isInter) for(k in 1:nIntra) C[,,k] <- M
      if(!isInter) for(j in 1:max(1,nInter)) C[,j,] <- M
      M <- C
      # Get rid of third dimension if nIntra=1
      if(dim(M)[3]==1) M <- M[,,1,drop=TRUE]
    }

    if(!silent) cat(".")
    return(M)
  }

  # Create list container for the draws
  drawsList <- list()

  # ############################################ #
  ### Load inter draws from global environment ###
  # ############################################ #
  if(testEUsr){
    if(!silent) cat("Reshaping inter-individual draws ")
    draws <- get(d$interDrawsType, envir=globalenv())

    # Validate input
    if( !is.list(draws) ) stop("Draws provided by user must be contained in a list.")
    for(i in draws ){
      if(!is.matrix(i)) stop("At least one element of ", d$interDrawsType, " is not a matrix.")
      if(nrow(i)!=nIndiv) stop("At least one element of ", d$interDrawsType, " does not have as many rows as individuals.")
      if(ncol(i)!=nInter) stop("At least one element of ", d$interDrawsType, " does not have ", nInter, " columns.")
    }
    
    # Turn uniform to standard normal if applicable
    for(i in d$interNormDraws) draws[[i]] <- stats::qnorm(draws[[i]])
    
    # Reshape draws into correct dimensionality
    drawsList <- c( drawsList, lapply(draws, expand, isInter=TRUE) )
    if(!silent) cat(" Done.\n")
  }

  # ############################################ #
  ### Load intra draws from global environment ###
  # ############################################ #
  if(testAUsr){
    if(!silent) cat("Reshaping intra-individual draws ")
    draws <- get(d$intraDrawsType, envir=globalenv())

    # Validate input
    if( !is.list(draws) ) stop("Draws provided by user must be contained in a list.")
    for(i in draws ){
      if(!is.matrix(i)) stop("At least one element of ", d$interDrawsType, " is not a matrix.")
      if(nrow(i)!=nObs) stop("At least one element of ", d$interDrawsType, " does not have as many rows as observations.")
      if(ncol(i)!=nIntra) stop("At least one element of ", d$intraDrawsType, " does not have ", nIntra, " columns.")
    }
    
    # Turn uniform to standard normal if applicable
    for(i in d$intraNormDraws) draws[[i]] <- stats::qnorm(draws[[i]])

    # Reshape draws into correct dimensionality
    drawsList <- c( drawsList, lapply(draws, expand, isInter=FALSE) )
    if(!silent) cat(" Done.\n")
  }

  # ######################################## #
  #### Generate inter draws               ####
  # ######################################## #
  if(all(testEGen)){
    if(!silent) cat("Generating inter-individual draws ")

    # Generate draws in a matrix. One column per dimension (i.e. random component)
    # and as many rows as nDraws*nIndiv
    if(d$interDrawsType=='halton') draws <- randtoolbox::halton(nInter*nIndiv, dimInter)
    if(d$interDrawsType=='mlhs') draws <- apollo_mlhs(nInter, dimInter, nIndiv)
    if(d$interDrawsType=='pmc') draws <- matrix(stats::runif(nIndiv*nInter*dimInter),
                                                            nrow=nIndiv*nInter, ncol=dimInter,
                                                            byrow=TRUE)
    if(d$interDrawsType=='sobol') draws <- randtoolbox::sobol(nInter*nIndiv, dimInter, scrambling=0)                                                            
    if(d$interDrawsType=='sobolowen') draws <- randtoolbox::sobol(nInter*nIndiv, dimInter, scrambling=1)
    if(d$interDrawsType=='sobolfauretezuka') draws <- randtoolbox::sobol(nInter*nIndiv, dimInter, scrambling=2)
    if(d$interDrawsType=='sobolowenfauretezuka') draws <- randtoolbox::sobol(nInter*nIndiv, dimInter, scrambling=3)
    
    draws <- as.matrix(draws)

    # Turn draws into a list of dimInter matrices of dimensions nIndiv x nInter.
    draws <- split(draws, rep(1:ncol(draws), each = nrow(draws)))
    draws <- lapply(draws, matrix, nrow=nIndiv, ncol=nInter, byrow=TRUE)
    names(draws) <- c(d$interUnifDraws, d$interNormDraws)

    # Turn uniform to standard normal if applicable
    for(i in d$interNormDraws) draws[[i]] <- stats::qnorm(draws[[i]])

    # Reshape draws into correct dimensionality
    drawsList <- c( drawsList, lapply(draws, expand, isInter=TRUE) )
    if(!silent) cat(" Done\n")
  }

  # ######################################## #
  #### Generate intra draws               ####
  # ######################################## #
  if(all(testAGen)){
    if(!silent) cat("Generating intra-individual draws ")

    # Generate draws in a matrix. One column per dimension (i.e. random component)
    # and as many rows as observations (nObs)
    if(d$intraDrawsType=='halton'){
      # Skip dimensions already used for inter
      haltonInter <- dimInter>0 && d$interDrawsType=='halton'
      draws <- randtoolbox::halton(nIntra*nObs, dimInter*haltonInter + dimIntra)
      if(haltonInter) draws <- draws[,(dimInter+1):(dimInter+dimIntra)]
    }
    if(d$intraDrawsType=='mlhs') draws <- apollo_mlhs(nIntra, dimIntra, nObs)
    if(d$intraDrawsType=='pmc') draws <- matrix(stats::runif(nIntra*nObs*dimIntra),
                                                 nrow=nIntra*nObs, ncol=dimIntra)
    sobolNames <- c('sobol', 'sobolowen', 'sobolfauretezuka', 'sobolowenfauretezuka')
    sobolInter <- dimInter>0 && (d$interDrawsType %in% sobolNames)
    if(d$intraDrawsType %in% sobolNames){
      sobolScrambling <- c(sobol=0, sobolOwen=1, sobolFaureTezuka=2, sobolOwenFaureTezuka=3)
      draws <- randtoolbox::sobol(nIntra*nObs, dimInter*sobolInter + dimIntra, 
                                  scrambling=sobolScrambling[d$intraDrawsType])
      if(sobolInter) draws <- draws[,(dimInter+1):(dimInter+dimIntra)]
    }
                                                
    draws <- as.matrix(draws)

    # Turn draws into a list of dimInter matrices of dimensions nIndiv x nInter.
    draws <- split(draws, rep(1:ncol(draws), each = nrow(draws)))
    draws <- lapply(draws, matrix, nrow=nObs, ncol=nIntra, byrow=TRUE)
    names(draws) <- c(d$intraUnifDraws, d$intraNormDraws)

    # Turn uniform to standard normal if applicable
    for(i in d$intraNormDraws) draws[[i]] <- stats::qnorm(draws[[i]])

    # Reshape draws into correct dimensionality
    drawsList <- c( drawsList, lapply(draws, expand, isInter=FALSE) )
    if(!silent) cat(" Done\n")
  }


  # ################################## #
  #### Returning draws              ####
  # ################################## #

  if(!panelData & dimInter>0) warning('Inter-person draws are being used without a panel structure.')

  return(drawsList)
}