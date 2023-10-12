#' Validates the \code{apollo_HB} list of parameters
#'
#' Validates the \code{apollo_HB} list of parameters and sets default values for the omitted ones.
#'
#' This function is only necessary when using bayesian estimation.
#' @param apollo_HB List. Contains options for Bayesian estimation. See \code{?RSGHB::doHB} for details.
#'                   Parameters \code{modelname}, \code{gVarNamesFixed}, \code{gVarNamesNormal},
#'                   \code{gDIST}, \code{svN} and \code{FC} are automatically set based on the
#'                   other arguments of this function.
#'                   Other settings to include are the following.
#'                   \itemize{
#'                     \item \strong{\code{constraintNorm}}: Character vector. Constraints for \emph{random} coefficients 
#'                                                          in bayesian estimation. Constraints can be written as 
#'                                                          "b1>b2", "b1<b2", "b1>0", or "b1<0".
#'                     \item \strong{\code{fixedA}}: Named numeric vector. Contains the names and fixed mean values of 
#'                                                  random parameters. For example, c(b1=0) fixes the mean of b1 to zero.
#'                     \item \strong{\code{fixedD}}: Named numeric vector. Contains the names and fixed variance of 
#'                                                  random parameters. For example, c(b1=1) fixes the variance of b1 to zero.
#'                     \item \strong{\code{gFULLCV:}} Boolean. Whether the full variance-covariance structure should be used for random parameters (TRUE by default).
#'                     \item \strong{\code{gNCREP}}: Numeric. Number of burn-in iterations to use prior to convergence (default=10^5).
#'                     \item \strong{\code{gNEREP}}: Numeric. Number of iterations to keep for averaging after convergence has been reached (default=10^5).
#'                     \item \strong{\code{gINFOSKIP}}: Numeric. Number of iterations between printing/plotting information about the iteration process (default=250).
#'                     \item \strong{\code{hbDist}}: \emph{Mandatory} setting. A named character vector determining
#'                                                  the distribution of each parameter to be estimated. Possible 
#'                                                  values are as follows.
#'                                                  \itemize{
#'                                                    \item \strong{\code{"CN+"}}: Positive censored normal.
#'                                                    \item \strong{\code{"CN-"}}: Negative censored normal.
#'                                                    \item \strong{\code{"JSB"}}: Johnson SB.
#'                                                    \item \strong{\code{"LN+"}}: Positive log-normal.
#'                                                    \item \strong{\code{"LN-"}}: Negative log-normal.
#'                                                    \item \strong{\code{"N"}}: Normal.
#'                                                    \item \strong{\code{"NR"}}: Fixed (as in non-random) parameter.
#'                                                  }
#'                     \item \strong{\code{nodiagnostics}}: Boolean. Turn off pre-estimation diagnostics for RSGHB. Set to TRUE by default.
#'                   }
#' @param apollo_beta Named numeric vector. Names and values for parameters.
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation.
#'                    value is constant throughout estimation).
#' @param apollo_control List. Options controlling the running of the code. See \link{apollo_validateInputs}.
#' @param silent Boolean. TRUE to keep the function from printing to the console. Default is FALSE.
#' @return Validated apollo_HB
#' @export
apollo_validateHBControl=function(apollo_HB, apollo_beta, apollo_fixed, apollo_control, silent=FALSE){
  
  if(is.null(apollo_HB$nodiagnostics)) apollo_HB$nodiagnostics=TRUE
  
  # Validate inputs
  if(!("hbDist" %in% names(apollo_HB))) stop("SYNTAX ISSUE - No hbDist element in apollo_HB.")
  
  if(any(apollo_HB$hbDist=="DNE")) apollo_HB$hbDist[which(apollo_HB$hbDist=="DNE")]="NR"
  
  hbDist <- apollo_HB$hbDist
  if(length(apollo_beta)!=length(hbDist)) stop("SYNTAX ISSUE - Argument hbDist has different length than apollo_beta.")
  if(!all(names(apollo_beta)==names(hbDist))) stop("SYNTAX ISSUE - Argument hbDist has differently named elements from apollo_beta (or in a different order).")
  
  hbDist_nonest = hbDist[(names(apollo_beta) %in% apollo_fixed)]
  if(any(!(hbDist_nonest%in%c("F","NR","DNE")))) stop("SYNTAX ISSUE - Only non-random parameters should be included in apollo_fixed for HB estimation, other constraints need to be accommodated in fixedA and fixedD!")
  
  # Translate rnd param distributions into RSGHB coding
  map <- c("NR"=0, "F"=0, "N"=1, "LN+"=2, "LN-"=3, "CN+"=4, "CN-"=5, "JSB"=6)
  hbDist_est = hbDist[!(names(apollo_beta) %in% apollo_fixed)]
  hbDist_est = stats::setNames(map[hbDist_est], names(hbDist_est)) # distribution in RSGHB coding
  theta_est  = apollo_beta[!(names(apollo_beta) %in% apollo_fixed)] # starting values
  
  
  gVarNamesFixed  = names(theta_est[hbDist_est==0])
  FC              = theta_est[hbDist_est==0]
  gVarNamesNormal = names(theta_est[hbDist_est>0])
  gDIST           = hbDist_est[hbDist_est>0]
  svN             = theta_est[hbDist_est>0]
  
  apollo_HB$modelname       = apollo_control$modelName
  apollo_HB$gVarNamesFixed  = gVarNamesFixed
  apollo_HB$gVarNamesNormal = gVarNamesNormal
  apollo_HB$gDIST           = gDIST
  apollo_HB$svN             = svN
  apollo_HB$FC              = FC
  
  if(is.null(apollo_HB$gFULLCV)) apollo_HB$gFULLCV=TRUE
  
  if(length(apollo_HB$gVarNamesNormal)==0) {
    apollo_HB$gVarNamesNormal = NULL
    apollo_HB$hIW             = FALSE
  }
  if(length(apollo_HB$gVarNamesFixed)==0) apollo_HB$gVarNamesFixed=NULL
  if(length(apollo_HB$gDIST)==0) apollo_HB$gDIST=NULL
  if(length(apollo_HB$svN)==0) apollo_HB$svN=NULL
  if(length(apollo_HB$FC)==0) apollo_HB$FC=NULL
  
  ### Translates fixedA to RSGHB coding
  fA <- apollo_HB$fixedA
  if(!is.null(fA) && is.vector(fA)  && is.numeric(fA) && !is.null(names(fA))){
    if(!all(names(fA) %in% apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - Some names in fixedA do not match the names of random parameters.")
    apollo_HB$fixedA <- rep(NA, length(apollo_HB$gVarNamesNormal))
    for(i in 1:length(fA)){
      j <- which(apollo_HB$gVarNamesNormal==names(fA)[i])
      apollo_HB$fixedA[j] <- fA[i]
    }; rm(i,j)
  }; rm(fA)
  
  ### Translates fixedD to RSGHB coding
  fD <- apollo_HB$fixedD
  if(!is.null(fD) && is.vector(fD)  && is.numeric(fD) && !is.null(names(fD))){
    if(!all(names(fD) %in% apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - Some names in fixedD do not match the names of random parameters.")
    apollo_HB$fixedD <- rep(NA, length(apollo_HB$gVarNamesNormal))
    for(i in 1:length(fD)){
      j <- which(apollo_HB$gVarNamesNormal==names(fD)[i])
      apollo_HB$fixedD[j] <- fD[i]
    }; rm(i,j)
  }; rm(fD)
  
  ### Additional checks
  if(!is.null(apollo_HB$fixedA) && length(apollo_HB$fixedA)!=length(apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - fixedA has a different length than gVarNamesNormal inside apollo_HB")
  if(!is.null(apollo_HB$fixedD) && length(apollo_HB$fixedD)!=length(apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - fixedD has a different length than gVarNamesNormal inside apollo_HB")
  if(any(hbDist_est==6)){
    if(is.null(apollo_HB$gMINCOEF) | is.null(apollo_HB$gMAXCOEF)) stop("SYNTAX ISSUE - JSB distribution in use, but gMINCOEF or gMAXCOEF not defined in apollo_HB")
    if(!is.null(apollo_HB$gMINCOEF) && length(apollo_HB$gMINCOEF)!=length(apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - gMINCOEF has a different length than gVarNamesNormal inside apollo_HB")
    if(!is.null(apollo_HB$gMAXCOEF) && length(apollo_HB$gMAXCOEF)!=length(apollo_HB$gVarNamesNormal)) stop("SYNTAX ISSUE - gMINCOEF has a different length than gVarNamesNormal inside apollo_HB")
  }
  
  ### Process constraints.
  if(!is.null(apollo_HB$constraintsNorm)){
    # If in old format, e.g.: list( c(1,1,2), c(3,2,5) )
    if(is.list(apollo_HB$constraintsNorm) && all(sapply(apollo_HB$constraintsNorm, is.numeric))){
      # Print constraints
      #nam <- names(theta_est)
      nam <- gVarNamesNormal
      if(!silent) apollo_print("Constraints:")
      for(i in 1:length(apollo_HB$constraintsNorm)){
        t0 <- apollo_HB$constraintsNorm[[i]]
        #test <- 0 < t0[1] & t0[1]<length(theta_est)
        test <- 0 < t0[1] & t0[1]<=length(nam)
        test <- test & t0[2] %in% 1:2
        #test <- test & ( 0 <= t0[3] & t0[3]<length(theta_est) )
        test <- test & ( 0 <= t0[3] & t0[3]<=length(nam) )
        if(!test) stop(paste0("SYNTAX ISSUE - Incorrect format for constraintsNorm element c(", t0, ")."))
        if(!silent) apollo_print(paste0(nam[t0[1]], ifelse(t0[2]==1, " < ", " > "), ifelse(t0[3]==0, 0, nam[t0[3]]), "\n", sep=""))
      }
    } else { # If in new format, e.g: c("b1<b2","b3>b5")
      # Validate input
      if(length(apollo_HB$constraintsNorm)<1) stop("SYNTAX ISSUE - constraintsNorm, if included, should be at least one element long")
      if(!is.character(apollo_HB$constraintsNorm)) stop("SYNTAX ISSUE - constraintsNorm, if included, should be a character vector")
      if(!all(grepl("[<>]", apollo_HB$constraintsNorm))) stop("SYNTAX ISSUE - constraintsNorm, if included, must contain constraints written as 'param1>param2', 'param1<param2', 'param1>0' or 'param1<0'")
      if(any(grepl("=", apollo_HB$constraintsNorm))) stop("SYNTAX ISSUE - constraintsNorm does not support equality (=), greater or equal (>=), or less or equal (<=) constraints")
      if(length(apollo_fixed)>0) for(f in apollo_fixed) if(any(grepl(f,apollo_HB$constraintsNorm))) stop("SYNTAX ISSUE - constraintsNorm, if included, should not include any fixed parameter")
      # Translate constraints
      translated <- list()
      apollo_HB$constraintsNorm <- gsub(" ", "", apollo_HB$constraintsNorm) # remove white spaces
      for(i in 1:length(apollo_HB$constraintsNorm)){
        con<- apollo_HB$constraintsNorm[i]
        t0 <- gregexpr("[<>]", con)[[1]][1]
        t1 <- substr(con, 1, t0-1)
        if(t1=="0") stop(paste0("SYNTAX ISSUE - Constraint ", con, " should have zero in the right hand side"))
        t2 <- substr(con, t0, t0)
        t3 <- substr(con, t0+1, nchar(con))
        t0 <- c(0,0,0)
        #if(!all(c(t1,t3) %in% c(gVarNamesFixed, gVarNamesNormal, "0"))) stop(paste0("SYNTAX ISSUE - Constraint ", con, " includes invalid variable names"))
        if(!all(c(t1,t3) %in% c(gVarNamesNormal, "0"))) stop(paste0("SYNTAX ISSUE - Constraint ", con, " includes invalid variable names (only random params allowed)"))
        #t0[1] <- which(names(theta_est)==t1)
        t0[1] <- which(gVarNamesNormal==t1)
        t0[2] <- ifelse(t2=="<", 1, 2)
        #t0[3] <- ifelse(t3=="0", 0, which(names(theta_est)==t3))
        t0[3] <- ifelse(t3=="0", 0, which(gVarNamesNormal==t3))
        translated[[i]] <- t0
      }
      names(translated) <- NULL
      apollo_HB$constraintsNorm <- translated
    }
  }
  
  if(!silent) apollo_print("All checks on apollo_HB completed.")
  return(apollo_HB)
}
