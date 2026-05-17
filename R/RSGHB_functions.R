# Functions adapted from RSGHB 1.2.2
# Original authors: Jeff Dumont, Jeff Keller, Chase Carpenter
# Original license: GPL-3
# Source: https://github.com/RSGInc/RSGHB
#' @importFrom graphics axis par plot points segments
#' @importFrom grDevices dev.flush graphics.off rainbow
#' @importFrom stats aggregate complete.cases density rnorm runif var
#' @importFrom utils write.table
#' @importFrom MCMCpack rinvgamma riwish
NULL
checkModel = function(nodiagnostics = FALSE, verbose = TRUE, env = parent.frame()) {
  
  passChecks = TRUE
  if (is.null(env$gDIST) & env$gNIV > 0) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Variable - gDIST - is undefined.\n")
  }
  if (is.null(env$gNCREP)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Variable - gNCREP - is undefined.\n")
  }
  if (is.null(env$gNEREP)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Variable - gNEREP - is undefined.\n")
  }
  if (is.null(env$gNSKIP)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Variable - gNSKIP - is undefined.\n")
  }
  if (is.null(env$gINFOSKIP)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Variable - gINFOSKIP - is undefined.\n")
  }
  if (is.null(env$likelihood)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: The likelihood function is undefined.\n")
  }
  if (env$gNIV + env$gFIV == 0) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Please specify at least one coefficient to be estimated in either in gVarNamesNormal or gVarNamesFixed.\n")
  }
  if (length(env$gDIST) != env$gNIV) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: The number of distributions specified in gDist doesn't equal the number of random coefficients in the model.\n")
  }
  if (env$gNIV > 0) {
    for (d in env$gDIST) {
      if (d < 1 | d > length(env$distNames)) {
        passChecks <- FALSE
        cat("\n********FATAL ERROR: The specified distributions ", 
            d, " in gDist do not exist\n")
      }
    }
  }
  if (env$gNIV != length(env$svN)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: There are too many/not enough starting values for the random coefficients. Check your sVN vector.\n")
  }
  if (env$gFIV != length(env$FC)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: There are too many/not enough starting values for the fixed coefficients. Check your FC vector.\n")
  }
  if (is.null(env$choicedata$ID)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: Expecting to find a respondent identifier column called - ID - in your dataset. None found.\n")
  }
  if (sum(sort(env$choicedata$ID) == env$choicedata$ID) != length(env$choicedata$ID)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: The choice data is not sorted by ID.\n")
  }
  if ((!is.null(env$fixedA)) & length(env$fixedA) != length(env$gVarNamesNormal)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: The fixedA vector is not of the same length as the gVarNamesNormal vector.\n")
  }
  if ((!is.null(env$fixedD)) & length(env$fixedD) != length(env$gVarNamesNormal)) {
    passChecks <- FALSE
    cat("\n********FATAL ERROR: The fixedD vector is not of the same length as the gVarNamesNormal vector.\n")
  }
  if (passChecks) prepareModel(env)
  if (passChecks & verbose) {
    cat(rep("\n", 5))
    cat("Diagnostic checks passed. Please review before proceeding\n")
    diagnostics <- data.frame(` ` = c("Number of Individuals:", 
                                      "Number of Observations:", "Custom Prior Matrix Used:", 
                                      "Prior variance:", "Target Acceptance (Fixed):", 
                                      "Target Acceptance (Normal):", "Degrees of Freedom:", 
                                      "Avg. Number of Observations per Individual:", "Initial Log-Likelihood:"), 
                              ` ` = as.character(rep(NA, 9)), check.names = FALSE, stringsAsFactors = FALSE)
    diagnostics[1, 2] <- env$gNP
    diagnostics[2, 2] <- env$gNOBS
    if (env$useCustomPVMatrix) {
      diagnostics[3, 2] <- TRUE
    } else {
      diagnostics[4, 2] <- signif(env$priorVariance, env$gSIGDIG)
    }
    if (env$gFIV > 0) 
      diagnostics[5, 2] <- signif(env$targetAcceptanceFixed, env$gSIGDIG)
    if (env$gNIV > 0) 
      diagnostics[6, 2] <- signif(env$targetAcceptanceNormal, env$gSIGDIG)
    diagnostics[7, 2] <- signif(env$degreesOfFreedom, env$gSIGDIG)
    diagnostics[8, 2] <- signif(env$gNOBS/env$gNP, env$gSIGDIG)
    diagnostics[9, 2] <- signif(sum(log(env$likelihood(env$FC, env$B, env))), env$gSIGDIG)
    cat("-----------------------------------------------------------\n")
    print(diagnostics[complete.cases(diagnostics), , drop = FALSE], row.names = FALSE)
    cat("\n-----------------------------------------------------------\n\n")
    if (env$gFIV > 0) {
      print(data.frame(`Fixed Parameters` = env$gVarNamesFixed, Start = env$FC, check.names = FALSE), row.names = FALSE)
      cat("\n-----------------------------------------------------------\n\n")
    }
    if (env$gNIV > 0) {
      print(data.frame(`Random Parameters` = env$gVarNamesNormal, 
                       Start = env$svN, Dist. = env$distNames[env$gDIST], 
                       check.names = FALSE), row.names = FALSE)
      cat("\n-----------------------------------------------------------\n\n")
    }
    if (!is.null(env$constraintsNorm)) {
      cat("Constraints applied to random parameters:\n")
      diagconstraints <- data.frame(` ` = unlist(lapply(env$constraintsNorm, 
                                                        FUN = `[`, 1)), ` ` = unlist(lapply(env$constraintsNorm, 
                                                                                            FUN = `[`, 2)), ` ` = unlist(lapply(env$constraintsNorm, 
                                                                                                                                FUN = `[`, 3)), check.names = FALSE)
      diagconstraints[, 1] <- env$gVarNamesNormal[diagconstraints[, 1]]
      diagconstraints[, 2] <- env$constraintLabels[diagconstraints[, 2]]
      diagconstraints[diagconstraints[, 3] != 0, 3] <- env$gVarNamesNormal[diagconstraints[diagconstraints[, 3] != 0, 3]]
      print(diagconstraints, row.names = FALSE, right = FALSE)
      cat("\n-----------------------------------------------------------")
    }
    if (!is.null(env$Choice)) {
      cat("\n", "Choice Matrix", "\n")
      choiceMatrix <- cbind(table(env$Choice), round(prop.table(table(env$Choice)), 2))
      dimnames(choiceMatrix)[[2]] <- c("Count", "%")
      print(choiceMatrix)
    }
    cat("\n\n\n")
    if (!nodiagnostics & verbose) {
      rl <- readline("Estimate Model? (Y/N): ")
      if (rl != "Y" & rl != "y") passChecks <- FALSE
    }
  }
  return(passChecks)
}

doHB = function(likelihood_user, choicedata, control = list()) {
  likelihood <- function(fc, b, env) {
    if (env$gNIV > 0) {
      gIDS <- env$gIDS
      C <- trans(b, env)
      if (env$gNIV > 1) 
        C <- C[gIDS, ]
      if (env$gNIV == 1) 
        C <- matrix(C[gIDS], ncol = env$gNIV)
    }
    p <- likelihood_user(fc, C)
    p <- replace(p, is.na(p), 9.88131291682493e-324)
    p0 <- rep(1, env$gNP)
    ## p0 <- apollo_aggregation(as.double(env$gIDS), as.double(env$gNOBS), 
    ##          as.double(env$gNP), as.double(p), as.double(1:env$gNP), 
    ##          as.double(p0))[[6]]
    ## avoiding C code
    if (any(p == 0)) {
      p0 <- as.numeric(tapply(p, env$gIDS, prod))
    } else {
      p0 <- exp(as.numeric(rowsum(log(p), env$gIDS, reorder = FALSE)))
    }
    p0 <- replace(p0, p0 < 9.88131291682493e-324, 9.88131291682493e-324)
    return(p0)
  }
  if (is.null(control[["modelname"]])) {
    modelname <- "HBModel"
  } else {
    modelname <- control[["modelname"]]
  }
  if (is.null(control[["gVarNamesNormal"]])) {
    gVarNamesNormal <- c()
  } else {
    gVarNamesNormal <- control[["gVarNamesNormal"]]
  }
  if (is.null(control[["gVarNamesFixed"]])) {
    gVarNamesFixed <- c()
  } else {
    gVarNamesFixed <- control[["gVarNamesFixed"]]
  }
  if (is.null(control[["gDIST"]])) {
    gDIST <- rep(1, length(gVarNamesNormal))
  } else {
    gDIST <- control[["gDIST"]]
  }
  if (is.null(control[["FC"]])) {
    FC <- rep(0, length(gVarNamesFixed))
  } else {
    FC <- control[["FC"]]
  }
  if (is.null(control[["svN"]])) {
    svN <- rep(0, length(gVarNamesNormal))
  } else {
    svN <- control[["svN"]]
  }
  if (is.null(control[["gNCREP"]])) {
    gNCREP <- 1e+05
  } else {
    gNCREP <- control[["gNCREP"]]
  }
  if (is.null(control[["gNEREP"]])) {
    gNEREP <- 1e+05
  } else {
    gNEREP <- control[["gNEREP"]]
  }
  if (is.null(control[["gNSKIP"]])) {
    gNSKIP <- 1
  } else {
    gNSKIP <- control[["gNSKIP"]]
  }
  if (is.null(control[["gINFOSKIP"]])) {
    gINFOSKIP <- 250
  } else {
    gINFOSKIP <- control[["gINFOSKIP"]]
  }
  if (is.null(control[["constraintsNorm"]])) {
    constraintsNorm <- NULL
  } else {
    constraintsNorm <- control[["constraintsNorm"]]
  }
  if (is.null(control[["fixedA"]])) {
    fixedA <- NULL
  } else {
    fixedA <- control[["fixedA"]]
  }
  if (is.null(control[["fixedD"]])) {
    fixedD <- NULL
  } else {
    fixedD <- control[["fixedD"]]
  }
  if (is.null(control[["nodiagnostics"]])) {
    nodiagnostics <- FALSE
  } else {
    nodiagnostics <- control[["nodiagnostics"]]
  }
  if (is.null(control[["gSIGDIG"]])) {
    gSIGDIG <- 10
  } else {
    gSIGDIG <- control[["gSIGDIG"]]
  }
  if (is.null(control[["priorVariance"]])) {
    priorVariance <- 2
  } else {
    priorVariance <- control[["priorVariance"]]
  }
  if (is.null(control[["degreesOfFreedom"]])) {
    degreesOfFreedom = 5
  } else {
    degreesOfFreedom <- control[["degreesOfFreedom"]]
  }
  ### CMC: additional checks for using hIW
  if (is.null(control[["hIW"]])|is.null(control[["gVarNamesNormal"]])) {
    hIW <- FALSE
  } else {
    hIW <- control[["hIW"]]
  }
  if (is.null(control[["xi"]]) & hIW == TRUE) {
    xi <- 1
  } else {
    xi <- control[["xi"]]
  }
  ### CMC: additional check to ensure > 1 when using hIW
  if (degreesOfFreedom < 2 & hIW == TRUE){
    degreesOfFreedom=2
    cat("Degrees of freedom set to 2 - lowest possible value")
  }
  if (is.null(control[["rho"]])) {
    rho <- 0.1
  } else {
    rho <- control[["rho"]]
  }
  if (is.null(control[["rhoF"]])) {
    rhoF <- 1e-04
  } else {
    rhoF <- control[["rhoF"]]
  }
  if (is.null(control[["gFULLCV"]])) {
    gFULLCV <- TRUE
  } else {
    gFULLCV <- control[["gFULLCV"]]
  }
  if (is.null(control[["gStoreDraws"]])) {
    gStoreDraws <- FALSE
  } else {
    gStoreDraws <- control[["gStoreDraws"]]
  }
  if (is.null(control[["gSeed"]])) {
    gSeed <- 0
  } else {
    gSeed <- control[["gSeed"]]
  }
  if (is.null(control[["gMINCOEF"]])) {
    gMINCOEF <- 0
  } else {
    gMINCOEF <- control[["gMINCOEF"]]
  }
  if (is.null(control[["gMAXCOEF"]])) {
    gMAXCOEF <- 0
  } else {
    gMAXCOEF <- control[["gMAXCOEF"]]
  }
  if (is.null(control[["pvMatrix"]])) {
    useCustomPVMatrix <- FALSE
    pvMatrix <- NULL
  } else {
    useCustomPVMatrix <- TRUE
    pvMatrix <- control[["pvMatrix"]]
  }
  if (is.null(control[["targetAcceptanceNormal"]])) {
    targetAcceptanceNormal <- 0.3
  } else {
    targetAcceptanceNormal <- control[["targetAcceptanceNormal"]]
  }
  if (is.null(control[["targetAcceptanceFixed"]])) {
    targetAcceptanceFixed <- 0.3
  } else {
    targetAcceptanceFixed <- control[["targetAcceptanceFixed"]]
  }
  if (is.null(control[["writeModel"]])) {
    writeModel <- FALSE
  } else {
    writeModel <- control[["writeModel"]]
  }
  if (is.null(control[["verbose"]])) {
    verbose <- TRUE
  } else {
    verbose <- control[["verbose"]]
  }
  
  gNP <- 0
  gNOBS <- 0
  TIMES <- 0
  gIDS <- 0
  respIDs <- 0
  A <- 0
  B <- 0
  Dmat <- 0
  Choice <- 0
  gNIV <- length(gVarNamesNormal)
  gFIV <- length(gVarNamesFixed)
  if (is.null(pvMatrix) & gNIV > 0) {
    pvMatrix <- priorVariance * diag(gNIV)
  }
  if (!is.matrix(pvMatrix) & gNIV > 0) {
    stop("\npvMatrix is not a matrix. Make sure that your prior covariance matrix is ", gNIV, " by ", gNIV, ".")
  }
  if (!is.null(pvMatrix) & gNIV > 0) {
    if (nrow(pvMatrix) != gNIV | ncol(pvMatrix) != gNIV) {
      stop("\nThe prior covariance matrix is of the wrong size. Make sure that your prior covariance matrix is ", gNIV, " by ", gNIV, ".")
    }
  }
  rownames(pvMatrix) <- colnames(pvMatrix) <- gVarNamesNormal
  
  begintime <- Sys.time()
  starttime <- Sys.time()
  distNames <- c("N", "LN+", "LN-", "CN+", "CN-", "JSB")
  constraintLabels <- c("<", ">")
  acceptanceRatePerc <- 0
  acceptanceRateF <- 0
  acceptanceRateFPerc <- 0
  rhoFadj <- 1e-05
  if (checkModel(nodiagnostics = nodiagnostics, verbose = verbose)) {
    r <- 1
    ma <- matrix(0, nrow = gNIV, ncol = gNEREP)
    md <- array(0, dim = c(gNIV, gNIV, gNEREP), dimnames = list(NULL, NULL, 1:gNEREP))
    mb <- matrix(0, nrow = gNP, ncol = gNIV)
    mb.squared <- matrix(0, nrow = gNP, ncol = gNIV)
    mp <- matrix(0, nrow = gNP, ncol = gNEREP)
    mf <- matrix(0, nrow = gFIV, ncol = gNEREP)
    mc <- matrix(0, nrow = gNP, ncol = gNIV)
    mc.squared <- matrix(0, nrow = gNP, ncol = gNIV)
    ### CMC: added two output matrices in two lines below
    cmcLLout  <- list()
    cmcRLHout <- matrix(0, nrow = gNP, ncol = gNEREP)
    storedDraws <- list()
    results <- list(modelname = modelname, params.fixed = gVarNamesFixed, 
                    params.vary = gVarNamesNormal, distributions = distNames[gDIST], 
                    pv = pvMatrix, df = degreesOfFreedom, gSIGDIG = gSIGDIG, 
                    gNP = gNP, gNOBS = gNOBS, gNCREP = gNCREP, gNEREP = gNEREP, 
                    gSeed = gSeed, constraints = constraintsNorm, iter.detail = data.frame(Iteration = NA, 
                                                                                           `Log-Likelihood` = NA, RLH = NA, `Parameter RMS` = NA, 
                                                                                           `Avg. Variance` = NA, `Acceptance Rate (Fixed)` = NA, 
                                                                                           `Acceptance Rate (Normal)` = NA, check.names = FALSE),
                    ### CMC: included new outputs in the results list in two lines below
                    cmcLLout = cmcLLout,
                    cmcRLHout = cmcRLHout)
    hb(A, B, Dmat, FC)
    matTIMES = matrix(1/TIMES,nrow=gNP,ncol=gNEREP,byrow = FALSE)
    if (gNIV > 0) {
      #browser()
      mp = cmcLLout^(matTIMES)
      ma <- cbind(iteration = (gNCREP + 1):(gNCREP + gNEREP), t(ma))
      mcsd <- cbind(id = respIDs, sqrt((mc.squared - mc^2/gNEREP)/gNEREP))
      mc <- cbind(id = respIDs, RLH = rowMeans(mp), mc/gNEREP)
      mbsd <- cbind(id = respIDs, sqrt((mb.squared - mb^2/gNEREP)/gNEREP))
      mb <- cbind(id = respIDs, mb/gNEREP)
      colnames(mc) <- c("Respondent", "RLH", gVarNamesNormal)
      colnames(ma) <- c("iteration", gVarNamesNormal)
      colnames(mcsd) <- c("Respondent", gVarNamesNormal)
      colnames(mb) <- c("Respondent", gVarNamesNormal)
      colnames(mbsd) <- c("Respondent", gVarNamesNormal)
      results$A <- ma
      results$B <- mb
      results$Bsd <- mbsd
      results$C <- mc
      results$Csd <- mcsd
      results$D <- md
    }
    if (gFIV > 0) {
      mf <- cbind(iteration = (gNCREP + 1):(gNCREP + gNEREP), t(mf))
      colnames(mf) <- c("iteration", gVarNamesFixed)
      results$F <- mf
    }
    ### CMC: independently of whether we've used random or fixed, we now save LL and RLH, next two lines added
    results$cmcLLout <- cmcLLout
    results$cmcRLHout <- results$cmcLLout^(matTIMES)
    if (gStoreDraws) {
      results$Draws <- storedDraws
      names(results$Draws) <- respIDs
    }
    results$choices <- Choice
    results$p <- likelihood_user(
      fc = if (is.null(results[["F"]])) {
        NULL
      } else {
        colMeans(as.matrix(results[["F"]][, -1]))
      }, 
      b = if (is.null(results[["C"]])) {
        NULL
      } else {
        as.matrix(results[["C"]][gIDS, -c(1:2)])
      })
    results$ll0 <- sum(log(likelihood_user(fc = FC, b = matrix(svN, ncol = length(svN), nrow = length(gIDS), byrow = TRUE))))
    results$llf <- sum(log(results[["p"]]))
    results[["iter.detail"]] <- results[["iter.detail"]][-1,]
    if (verbose) {
      cat("Estimation complete.\n")
    }
    class(results) <- "RSGHB"
    if (writeModel) {
      if (verbose) {
        cat("Creating output files. Please be patient.\n")
      }
      writeRSGHBModel(results)
      if (verbose) {
        cat("Output files finished writing to working directory.\n")
      }
    }
  } else {
    results <- NULL
  }
  return(results)
}

hb = function(a, b, d, f, env = parent.frame()) {
  
  env$starttime <- Sys.time()
  
  if (env$gSeed == 0) {
    env$gSeed <- ceiling(runif(1) * 1e+06)
  }
  set.seed(env$gSeed, kind = "default", normal.kind = "default")
  
  ### CMC: set lambda 
  if(env$hIW){
    lambda = diag(env$gNIV)
    for(l in 1:env$gNIV){
      lambda[l,l] = 1/rinvgamma(1,shape=(1/2),scale=(1/(env$xi^2)))
    }
  }
  
  p <- env$likelihood(f, b, env)
  
  for (r in 1:env$gNCREP) {
    if (env$gNIV > 0) {
      out <- nextB(a, b, d, p, f, env)
      b <- out[[1]]
      env$acceptanceRatePerc <- out[[2]]
      p <- out[[3]]
      a <- nextA(b, d, env)
      if (!is.null(env$fixedA)) {
        for (rp in 1:env$gNIV) {
          if (!is.na(env$fixedA[rp])) {
            a[rp] <- env$fixedA[rp]
          }
        }
      }
      if (env$gFULLCV){ 
        ### CMC: different functions used depending on whether hIW is true
        if(env$hIW){
          d <- nextDhIW(a, b, env,lambda)
          lambda <- nextlambda(d,env)
        } else {
          d <- nextD(a, b, env)  
        }
      }
      if (!env$gFULLCV) {
        d <- nextDind(a, b, env)
      }
      if (!is.null(env$fixedD)) {
        for (rp in 1:env$gNIV) {
          if (!is.na(env$fixedD[rp])) {
            d[rp, rp] <- env$fixedD[rp]
          }
        }
      }
    }
    if (env$gFIV > 0) {
      out <- nextF(p, f, b, env)
      if (sum(out[[1]] == f) != env$gFIV) {
        env$acceptanceRateF <- env$acceptanceRateF + 1
      }
      if ((r%%100) == 0) {
        env$acceptanceRateFPerc <- env$acceptanceRateF/100
        if (env$acceptanceRateFPerc < env$targetAcceptanceFixed) {
          env$rhoF <- env$rhoF - env$rhoF/50
        }
        if (env$acceptanceRateFPerc > env$targetAcceptanceFixed) {
          env$rhoF <- env$rhoF + env$rhoF/50
        }
        env$acceptanceRateF <- 0
      }
      f <- out[[1]]
      p <- out[[2]]
    }
    if (r%%env$gINFOSKIP == 0 | r == 1) {
      progreport(r, p, a, b, d, f, env)
    }
    if (env$gStoreDraws) {
      for (i in 1:env$gNP) {
        env$storedDraws[[i]] <- matrix(0, env$gNEREP, env$gNIV)
        dimnames(env$storedDraws[[i]]) <- list(NULL, env$gVarNamesNormal)
      }
    }
  }
  #browser()
  n <- env$gNEREP * env$gNSKIP
  for (r in 1:n) {
    if (env$gNIV > 0) {
      out <- nextB(a, b, d, p, f, env)
      b <- out[[1]]
      env$acceptanceRatePerc <- out[[2]]
      p <- out[[3]]
      a <- nextA(b, d, env)
      if (!is.null(env$fixedA)) {
        for (rp in 1:env$gNIV) {
          if (!is.na(env$fixedA[rp])) {
            a[rp] <- env$fixedA[rp]
          }
        }
      }
      if (env$gFULLCV){ 
        ### CMC: different functions used depending on whether hIW is true
        if(env$hIW){
          d <- nextDhIW(a, b, env,lambda)
          lambda <- nextlambda(d,env)
        } else {
          d <- nextD(a, b, env)  
        }
      }
      if (!env$gFULLCV) {
        d <- nextDind(a, b, env)
      }
      if (!is.null(env$fixedD)) {
        for (rp in 1:env$gNIV) {
          if (!is.na(env$fixedD[rp])) {
            d[rp, rp] <- env$fixedD[rp]
          }
        }
      }
      if (r%%env$gNSKIP == 0) {
        C <- trans(b, env)
        env$ma[, r/env$gNSKIP] <- a
        env$md[, , r/env$gNSKIP] <- d
        #env$mp[, r/env$gNSKIP] <- p^(1/env$TIMES)
        env$mb <- env$mb + b
        env$mb.squared <- env$mb.squared + b^2
        env$mc <- env$mc + C
        env$mc.squared <- env$mc.squared + C^2
        
        if (env$gStoreDraws) {
          for (i in 1:env$gNP) {
            env$storedDraws[[i]][(r/env$gNSKIP), ] <- C[i, ]
          }
        }
      }
    }
    if (env$gFIV > 0) {
      
      out <- nextF(p, f, b, env)
      if (sum(out[[1]] == f) != env$gFIV) {
        env$acceptanceRateF <- env$acceptanceRateF + 1
        if ((r%%100) == 0) {
          env$acceptanceRateFPerc <- env$acceptanceRateF/100
          if (env$acceptanceRateFPerc < env$targetAcceptanceFixed) {
            env$rhoF <- env$rhoF - env$rhoF/50
          }
          if (env$acceptanceRateFPerc > env$targetAcceptanceFixed) {
            env$rhoF <- env$rhoF + env$rhoF/50
          }
          env$acceptanceRateF <- 0
        }
      }
      f <- out[[1]]
      p <- out[[2]]
      
      if (r %% env$gNSKIP == 0) {
        env$mf[, (r/env$gNSKIP)] <- f
      }
    }
    ### CMC: independently of whether we've used random or fixed, we now save LL and RLH, next four lines added
    if (r %% env$gNSKIP == 0) {
      env$cmcLLout[[(r / env$gNSKIP)]] <- p
      #env$cmcRLHout[, (r/env$gNSKIP)] <- p^(1/env$TIMES)
    }
    if (r %% env$gINFOSKIP == 0) {
      progreport(env$gNCREP + r, p, a, b, d, f, env)
    }
  }
  env$cmcLLout = do.call(cbind, env$cmcLLout)
  return(TRUE)
}


### plot.RSGHB <- function(x, ...) { # add column argument?
###   
###   # Store old graphical parameters for later
###   old.par <- par(no.readonly = TRUE)
###   
###   if (is.null(as.list(match.call())$type)) {
###     type <- "Log"
###   } else {
###     type <- as.list(match.call())$type
###   }
###   
###   # Plot means
###   if (type == "A" | type == "F") {
###     A <- x[[type]]
###     
###     if (is.null(A)) stop(paste0("model object does not contain component ", type))
###     
###     p <- ncol(A) - 1
###     
###     # Arrange plots in a roughly square grid
###     par(oma = c(0, 0, 2, 0)) # can the margins be tightened further?
###     if (p < 4) {
###       par(mfrow = c(p, 1))
###     } else {
###       r <- ceiling(sqrt(p))
###       if (r * (r - 1) >= p) {c <- r - 1} else {c <- r}               
###       par(mfrow = c(r, c))
###     }
###     
###     # Plot
###     for (i in 1:p) plot(x = A[, 1], y = A[, 1 + i], type = "l", xlab = "Iteration", ylab = "Estimate", main = colnames(A)[1 + i])
###     mtext("Markov Chains", outer = TRUE, cex = 1.5)          
###     
###     #      } else if (type == "B") {
###     #           
###     #      } else if (type == "C") {
###     
###   } else if (type == "Log") {
###     
###     logStats <- x[["iter.detail"]]
###     
###     # Get valid columns
###     stats <- names(logStats)[-1]
###     valid.stats <- c()
###     for (stat in stats) {
###       if (!is.null(logStats[, stat]) & !all(is.na(logStats[, stat]))) valid.stats <- c(valid.stats, stat)
###     }
###     
###     # Plot
###     par(mfrow = c(length(valid.stats), 1), mar = c(4.1, 4.1, 2.1, 2.1))
###     for (stat in valid.stats) {
###       plot(x = logStats[, "Iteration"], logStats[, stat], type = "l", xlab = "Iteration", ylab = stat)
###     }
###   } else {
###     
###     stop("Invalid 'type' argument")
###     
###   }
###   
###   # Restore old graphical parameters
###   par(old.par)
###   
### }
### 
### 
### print.RSGHB <- function(x,...) {
###   model = x
###   
###   cat("Model:", model[["modelname"]])
###   cat("\n\n")
###   cat("Individuals:", length(unique(model[["C"]][, "Respondent"])))
###   cat("\n")
###   cat("Iterations Kept:", nrow(model[["A"]]))
###   cat("\n\n")
###   
###   # summary fit statistics
###   # report statistics for posterior iterations
###   cat("Fit statistics\n")
###   posterior <- (model[["iter.detail"]]$Iteration > model[["gNCREP"]])
###   
###   # need to make this clear if it is the lower level model
###   if(!is.null(model[["A"]])){
###     #cat("Mean log-likelihood (upper-level):",mean(model$sLikelihood),"\n")     
###   }
###   
###   cat("Mean log-likelihood (lower-level):", mean(model[["iter.detail"]][posterior,"Log-Likelihood"]),"\n")
###   cat("Mean root likelihood:", mean(model[["iter.detail"]][posterior,"RLH"]),"\n")
###   #cat("Hit rate:","***NEEDS Predict method**** hit rate table\n")
###   cat("\n\n")
###   
###   #cat("Model comparisons statistics\n")
###   #cat("Bayes Factor:",mean(model$sLikelihood)/model[["ll0"]],"\n") # can be used in calculating bayes factor
###   #D_hat <- -2*mean(model$sLikelihood)
###   #p_d   <- D_hat # - D(mean)
###   #DIC   <- p_d + D_hat
###   
###   #cat("Deviance Information Criterion (not complete):",DIC,"\n")
###   cat("\n")
###   
###   # if has random parameters
###   if(!is.null(model[["A"]])){
###     posterior.means <- colMeans(model[["A"]][,-1])
###     posterior.stdev <- apply(model[["A"]][,-1],2,sd)
###     
###     cat("Random Parameters (Underlying Normals)\n")
###     cat("---------------------------------------------\n")
###     cat("                        95% Credible Regions \n")
###     print(
###       data.frame(
###         Estimate  = signif(posterior.means, 3),
###         `Std Dev` = signif(posterior.stdev, 3),
###         `Min`     = apply(model[["A"]][,-1],2,quantile,0.025),
###         `Max`     = apply(model[["A"]][,-1],2,quantile,0.975),
###         check.names = FALSE,
###         row.names = model[["params.vary"]]
###       )
###     )
###   }
###   
###   # if has fixed parameters
###   if(!is.null(model[["F"]])){
###     posterior.means <- colMeans(model[["F"]][,-1])
###     posterior.stdev <- apply(model[["F"]][,-1],2,sd)
###     
###     cat("Fixed Parameters\n")
###     cat("---------------------------------------------\n")
###     cat("                        95% Credible Regions \n")
###     print(
###       data.frame(
###         Estimate  = signif(posterior.means, 3),
###         `Std Dev` = signif(posterior.stdev, 3),
###         `Min`     = apply(model[["F"]][,-1],2,quantile,0.025),
###         `Max`     = apply(model[["F"]][,-1],2,quantile,0.975),
###         check.names = FALSE,
###         row.names = model[["params.fixed"]]
###       )
###     )
###   }     
###   
### }
### 
### summary.RSGHB <- function(object,...) {
###   model = object
###   print(model)
### }

nextA = function(b, d, env) {
  SQNP <- sqrt(env$gNP)
  
  # draws from multivariate normals can be taken using 
  # draws = mu + l * eta where l*l' = d (l is the cholesky decomposition)      
  
  # accounting for non-diffuse priors
  newDraw = colMeans(b) + t(chol(d)) %*% matrix(rnorm(env$gNIV), nrow = env$gNIV, ncol = 1)/SQNP 
  return(newDraw)
}

nextB = function(a, b, d, p, f, env) {
  
  constraintsNorm <- env$constraintsNorm
  
  # note the multiplication by rho. this is the "jumping distribution"
  bnew <- b + matrix(rnorm(env$gNP * env$gNIV), nrow = env$gNP, ncol = env$gNIV) %*% chol(d) * sqrt(env$rho)
  
  # apply constraints here
  # cycle through constraints till no constraints are violated
  if (!is.null(constraintsNorm)) {
    totalConstraintsViolated <- 999
    while (totalConstraintsViolated > 0) {
      totalConstraintsViolated <- 0
      for (i in 1:length(constraintsNorm)) {
        firstpar <- constraintsNorm[[i]][1]
        condition <- constraintsNorm[[i]][2]
        secondpar <- constraintsNorm[[i]][3]
        if (condition == 1) {
          if (secondpar == 0) {
            constraintsViolated <- bnew[, firstpar] > 0
            bnew[constraintsViolated, firstpar] <- 0
          }
          if (secondpar != 0) {
            constraintsViolated <- bnew[, firstpar] > bnew[, secondpar]
            bnew[constraintsViolated, c(firstpar, secondpar)] <- bnew[constraintsViolated, secondpar]
          }
        }
        if (condition == 2) {
          if (secondpar == 0) {
            constraintsViolated <- bnew[, firstpar] < 0
            bnew[constraintsViolated, firstpar] <- 0
          }
          if (secondpar != 0) {
            constraintsViolated <- bnew[, firstpar] < bnew[, secondpar]
            bnew[constraintsViolated, c(firstpar, secondpar)] <- bnew[constraintsViolated, secondpar]
          }
        }
        totalConstraintsViolated <- totalConstraintsViolated + sum(constraintsViolated)
      }
    }
  }
  bn <- bnew - matrix(t(a), nrow = env$gNP, ncol = env$gNIV, byrow = T)
  bo <- b - matrix(t(a), nrow = env$gNP, ncol = env$gNIV, byrow = T)
  pnew <- env$likelihood(f, bnew, env)
  
  # gives the relative density (or the probability of seeing the betas) given current estimates
  # of D and A	
  # the relative densities are important because without any reference to how the betas fit in with the current
  # estimates of A and D, you'd end up with 
  
  r.new <- (log(pnew) + -0.5 * (colSums(t(bn) * (solve(d) %*% t(bn)))) - log(p) - -0.5 * (colSums(t(bo) * (solve(d) %*% t(bo)))))
  
  # if r.new > 1 then we accept the new estimate of beta. if r.new < 1 then we accept the new estimate
  # with probability = r.new
  ind <- (r.new >= 0) + (r.new < 0) * (log(matrix(runif(env$gNP), nrow = env$gNP)) <= r.new)
  nind <- 1 * (ind == 0)
  i <- colSums(ind)/env$gNP
  
  
  # this is the acceptance rate. the target for this 0.3 (though Sawtooth allows for the user to specify this).
  if (i < env$targetAcceptanceNormal) {
    env$rho <- env$rho - env$rho/10
  }
  if (i > env$targetAcceptanceNormal) {
    env$rho <- env$rho + env$rho/10
  }
  
  # I've just converted it to matrix form to make the multiplication simpler.
  mind <- matrix(0, nrow = env$gNP, ncol = env$gNIV)
  mnind <- matrix(0, nrow = env$gNP, ncol = env$gNIV)
  mind[, 1:env$gNIV] <- ind
  mnind[, 1:env$gNIV] <- 1 * (ind == 0)
  
  return(list(mind * bnew + mnind * b, i, ind * pnew + nind * p))
}


nextD = function(a, b, env) {
  
  `a-b` <- matrix(a, nrow = env$gNP, ncol = env$gNIV, byrow = TRUE) - b
  
  H <- env$pvMatrix + t(`a-b`) %*% `a-b`
  
  Tmat <- t(chol(solve(H)))
  
  u <- matrix(rnorm(env$gNIV * (env$gNP + env$gNIV + env$degreesOfFreedom)), 
              nrow = env$gNIV, ncol = env$gNP + env$gNIV + env$degreesOfFreedom)
  
  S <- (Tmat %*% u) %*% t(Tmat %*% u)
  
  return(solve(S))
}

nextDhIW = function(a, b, env, lambda) {
  
  S = (b - t(matrix(a,nrow=env$gNIV,ncol=env$gNP)))
  
  Vnew = t(S)%*%S + 2*env$degreesOfFreedom*lambda
  
  Sigmanew = riwish(env$degreesOfFreedom+env$gNP+env$gNIV-1,Vnew)
  
  return(Sigmanew)
}

nextDind = function(a, b, env) {
  
  d <- matrix(0, env$gNIV, env$gNIV)
  b <- matrix(t(a), nrow = env$gNP, ncol = env$gNIV, byrow = TRUE) - b
  
  for (k in 1:env$gNIV) {
    t <- 1 + t(b[, k]) %*% b[, k]
    s <- sqrt(1/t)[1, 1] * matrix(rnorm(env$gNP + 1), nrow = env$gNP + 1, ncol = 1)
    d[k, k] <- solve(t(s) %*% s)
  }
  
  return(d)
}

nextF = function(p, f, b, env) {
  
  # applying an MH algo here  
  fnew <- f + sqrt(env$rhoF) * rnorm(env$gFIV, 0, 1)
  pnew <- env$likelihood(fnew, b, env)
  
  r <- sum(log(pnew) - log(p))
  ind <- (log(runif(1)) <= r)
  
  pnew <- pnew * ind + p * (1 - ind)
  fnew <- as.vector(t(fnew %*% t(ind)) + t(f %*% t(1 - ind)))
  
  return(list(fnew, pnew))
}

nextlambda = function(d,env){
  
  lambdanew = diag(env$gNIV)
  
  invSigma = chol2inv(chol(d))
  
  for(j in 1:env$gNIV){
    lambdanew[j,j] = 1/rinvgamma(1,shape=((env$degreesOfFreedom+env$gNIV)/2),scale=((1/(env$xi^2))+(env$degreesOfFreedom*invSigma[j,j])))
  }
  
  return(lambdanew)
}


plotC = function(object, columns = NULL) {
  
  if (is.null(object[["C"]])) {
    stop("No random parameters to plot.")
  }
  
  C <- object[["C"]]
  gNIV <- ncol(C) - 2
  numIDs <- nrow(C)
  graphics.off()
  par(ask = TRUE)
  
  if (is.null(columns)) {
    columns <- 1:gNIV + 2
  }
  for (i in columns) {
    plot(density(C[, i]), type = "l", main = paste(colnames(C)[i], 
                                                   "\n", sum(C[, i] >= 0), " (", round(sum(C[, i] >= 0)/length(C[, i]) * 100, 1), "%) >= 0", sep = ""), 
         xlab = "Utility", ylab = "Density")
    dev.flush()
  }
  par(ask = FALSE)
}

prepareModel = function(env) {
  
  env$gNP <- length(unique(env$choicedata$ID))
  env$Choice <- env$choicedata$Choice
  env$gNOBS <- dim(env$choicedata)[1]
  env$TIMES <- matrix(0, nrow = env$gNP, ncol = 1)
  env$TIMES[, 1] <- aggregate(env$choicedata$ID, by = list(env$choicedata$ID), length)[, 2]
  env$gIDS <- unlist(as.vector(mapply(rep, 1:env$gNP, env$TIMES)))
  env$respIDs <- unique(env$choicedata$ID)
  
  if (length(env$gVarNamesNormal) > 0) {
    env$A <- matrix(0, nrow = env$gNIV, ncol = 1)
    env$B <- matrix(0, nrow = env$gNP, ncol = env$gNIV)
    env$Dmat <- env$priorVariance * diag(env$gNIV)
    env$A[, 1] <- env$svN
    env$B <- 1 + env$B
    env$B <- env$B * matrix(t(env$A), nrow = env$gNP, ncol = env$gNIV, byrow = T)
  }
}

progreport = function(r, p, a, b, d, f, env) {
  
  if (env$gNIV > 0) {
    paramRMS <- sqrt(mean(apply(trans(b, env), 1, function(x) x^2)))
    avgVariance <- mean(apply(trans(b, env), 1, var))
  } else {
    paramRMS <- NA
    avgVariance <- NA
  }
  if (env$verbose) {
    cat(rep("\n", 5))
    cat("-----------------------------------------------------------\n")
    cat("Iteration: ", r, "\n", sep = "\t")
    cat("-----------------------------------------------------------\n")
    mstats <- data.frame(
      ` ` = c("RHO (Fixed):", "Acceptance Rate (Fixed):", 
              "RHO (Normal):", "Acceptance Rate (Normal):", "Parameter RMS:", 
              "Avg. Variance:", "Log-Likelihood:", "RLH:"), ` ` = as.character(NA), 
      check.names = FALSE, stringsAsFactors = FALSE)
    
    if (env$gFIV > 0) {
      mstats[1, 2] <- signif(env$rhoF, env$gSIGDIG)
      mstats[2, 2] <- signif(env$acceptanceRateFPerc, env$gSIGDIG)
    }
    if (env$gNIV > 0) {
      mstats[3, 2] <- signif(env$rho, env$gSIGDIG)
      mstats[4, 2] <- signif(env$acceptanceRatePerc, env$gSIGDIG)
      mstats[5, 2] <- signif(paramRMS, env$gSIGDIG)
      mstats[6, 2] <- signif(avgVariance, env$gSIGDIG)
    }
    mstats[7, 2] <- signif(sum(log(p)), env$gSIGDIG)
    mstats[8, 2] <- signif(mean(p^(1/env$TIMES)), env$gSIGDIG)
    print(mstats[complete.cases(mstats), , drop = FALSE], right = TRUE, row.names = FALSE)
    cat("\n-----------------------------------------------------------\n\n")
    if (env$gFIV > 0) {
      print(data.frame(`Fixed Parameters` = paste0(env$gVarNamesFixed, ":"), Estimate = signif(f, env$gSIGDIG), check.names = FALSE), row.names = FALSE)
      cat("\n-----------------------------------------------------------\n\n")
    }
    if (env$gNIV > 0) {
      print(data.frame(`Random Parameters` = paste0(env$gVarNamesNormal, ":"), Estimate = signif(a, env$gSIGDIG), check.names = FALSE), row.names = FALSE)
      cat("\n-----------------------------------------------------------\n")
    }
    if (r > 1) {
      tpi <- (Sys.time() - env$starttime)/env$gINFOSKIP
      env$starttime <- Sys.time()
      tleft <- (env$gNCREP + env$gNEREP * env$gNSKIP - r) * tpi
      units(tleft) <- "mins"
      cat("Time per iteration:", format(tpi, digits = 3))
      cat("\n")
      cat("Time to completion:", format(tleft, digits = 3))
      cat("\n")
    } else {
      cat("Time per iteration: Calculating...")
      cat("\n")
      cat("Time to completion: Calculating...")
      cat("\n")
    }
    cat("-----------------------------------------------------------\n")
    if (env$gNIV > 0 & env$gFIV > 0) {
      alphas <- c(r, a, f)
    } else if (env$gNIV > 0) {
      alphas <- c(r, a)
    } else if (env$gFIV > 0) {
      alphas <- c(r, f)
    }
    
    cr <- rainbow(length(alphas) - 1)
    if (r == 1) {
      xmax <- (env$gNCREP + env$gNEREP * env$gNSKIP) * 1.05
      plot(x = 0, y = 0, main = "Markov Chains", xlim = c(0, xmax), ylim = c(-5, 5), pch = 20, xlab = "Iterations", ylab = "Utility", axes = FALSE, col = "white", cex = 0.5)
      segments(env$gNCREP, -100, env$gNCREP, 100, col = "red", lty = 2, lwd = 2)
      segments(0, 0, env$gNCREP + env$gNEREP * env$gNSKIP, 0, col = "gray", lty = 1, lwd = 1)
      axis(1, at = seq(from = 0, to = env$gNCREP + env$gNEREP * env$gNSKIP, by = floor((env$gNCREP + env$gNEREP * env$gNSKIP)/10)))
      axis(2, at = -100:100)
    }
    for (i in 2:length(alphas)) {
      points(x = alphas[1], y = alphas[i], pch = 20, col = cr[i - 1], cex = 0.5)
    }
    Sys.sleep(0)
  }
  detail <- c(r, signif(sum(log(p)), env$gSIGDIG), signif(mean(p^(1/env$TIMES)), 
                                                          env$gSIGDIG), signif(paramRMS, env$gSIGDIG), signif(avgVariance, 
                                                                                                              env$gSIGDIG), 
              if (env$gFIV > 0) {
                signif(env$acceptanceRateFPerc, env$gSIGDIG)
              } else {
                NA
              }, 
              if (env$gNIV > 0) {
                signif(env$acceptanceRatePerc, env$gSIGDIG)
              } else {
                NA
              })
  
  env$results[["iter.detail"]] <- rbind(env$results[["iter.detail"]], detail)
}

trans = function(b, env) {
  C <- b
  gDIST <- env$gDIST
  gMAXCOEF <- env$gMAXCOEF
  gMINCOEF <- env$gMINCOEF
  for (k in 1:dim(b)[2]) {
    if (gDIST[k] == 2) {
      C[, k] <- exp(b[, k])
    }
    if (gDIST[k] == 3) {
      C[, k] <- -1 * exp(b[, k])
    }
    if (gDIST[k] == 4) {
      C[, k] <- b[, k] * (b[, k] >= 0)
    }
    if (gDIST[k] == 5) {
      C[, k] <- b[, k] * (b[, k] <= 0)
    }
    if (gDIST[k] == 6) {
      C[, k] <- exp(b[, k])/(1 + exp(b[, k]))
      C[, k] <- (C[, k] * (gMAXCOEF[k] - gMINCOEF[k])) + gMINCOEF[k]
    }
  }
  return(C)
}

writeRSGHBModel = function(object, writeDraws = FALSE, path = getwd()) {
  
  options(max.print=1000000)     
  
  if (!"RSGHB" %in% class(object)) {
    stop("'object' is not of class RSGHB")
  }
  orig.path <- getwd()
  setwd(path)
  gSIGDIG <- object[["gSIGDIG"]]
  modelname <- object[["modelname"]]
  orig <- modelname
  i <- 1
  while (any(file.exists(paste0(modelname, c(".log", "_A.csv", "_B.csv", "_Bsd.csv", "_C.csv", "_Csd.csv", "_D.csv", "_F.csv"))))) {
    modelname <- paste0(orig, "~", i)
    i <- i + 1
  }
  if (modelname != orig) {
    warning(paste0("Model files associated with '", orig, "' already exist. Writing results as '", modelname, "'"))
  }
  orig.options <- options()
  options(width = 1000)
  sink(paste0(object[["modelname"]], ".log"))
  cat("Model Name:", object[["modelname"]], "\n", sep = "\t")
  cat("Number of individuals:", object[["gNP"]], "\n", sep = "\t")
  cat("Number of observations:", object[["gNOBS"]], "\n", sep = "\t")
  cat("Number of preliminary iterations:", object[["gNCREP"]], "\n", sep = "\t")
  cat("Number of draws used per individual:", object[["gNEREP"]], "\n", sep = "\t")
  cat("Random Seed:", object[["gSeed"]], "\n", sep = "\t")
  cat("Total iterations:", object[["gNCREP"]] + object[["gNEREP"]], "\n", sep = "\t")
  if (!is.null(object[["df"]])) {
    cat("Degrees of Freedom:", object[["df"]], "\n", sep = "\t")
  }
  cat("Number of parameters:", length(object$params.vary) + length(object$params.fixed), "\n\n", sep = "\t")
  if (length(object[["params.fixed"]]) > 0) {
    cat("Fixed parameters estimated:\n")
    cat(paste0(object[["params.fixed"]], "\n", collapse = ""))
  }
  cat("\n")
  if (length(object[["params.vary"]]) > 0) {
    cat("Random parameters estimated (Distribution):", "\n")
    cat(paste0(paste(object[["params.vary"]], "(", object[["distributions"]], ")"), "\n"), collapse = "", sep = "")
  }
  cat("\n")
  if (!is.null(object[["constraints"]])) {
    cond <- c("<", ">")
    cat("Constraints applied to random parameters (param1 - inequality - param2):\n")
    for (i in 1:length(object[["constraints"]])) {
      if (object[["constraints"]][[i]][3] == 0) {
        cat(object[["params.vary"]][object[["constraints"]][[i]][1]], cond[object[["constraints"]][[i]][2]],0, "\n")
      }
      if (object[["constraints"]][[i]][3] != 0) {
        cat(object[["params.vary"]][object[["constraints"]][[i]][1]], cond[object[["constraints"]][[i]][2]], object[["params.vary"]][object[["constraints"]][[i]][3]], "\n")
      }
    }
  }
  cat("\n")
  if (!is.null(object[["pv"]])) {
    cat("Prior Variance-Covariance Matrix:\n")
    print(object[["pv"]])
  }
  cat("\n-----------------------------------------------------------\n\n")
  print(object[["iter.detail"]], row.names = FALSE)
  sink()
  options(orig.options)
  if (!is.null(object[["pv"]])) {
    write.table(object[["pv"]], paste0(modelname, "_pvMatrix.csv"), sep = ",", col.names = NA)
  }
  if (!is.null(object[["A"]])) {
    write.table(signif(object[["A"]], gSIGDIG), paste0(modelname, "_A.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["B"]])) {
    write.table(signif(object[["B"]], gSIGDIG), paste0(modelname, "_B.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["Bsd"]])) {
    write.table(signif(object[["Bsd"]], gSIGDIG), paste0(modelname, "_Bsd.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["C"]])) {
    write.table(signif(object[["C"]], gSIGDIG), paste0(modelname, "_C.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["Csd"]])) {
    write.table(signif(object[["Csd"]], gSIGDIG), paste0(modelname, "_Csd.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["F"]])) {
    write.table(signif(object[["F"]], gSIGDIG), paste0(modelname, "_F.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["D"]])) {
    labelmatrix <- matrix(1:(length(object[["params.vary"]])^2), 
                          length(object[["params.vary"]]), length(object[["params.vary"]]))
    rownames(labelmatrix) <- colnames(labelmatrix) <- object[["params.vary"]]
    dlabels <- paste0(rownames(labelmatrix)[row(labelmatrix)[lower.tri(labelmatrix, diag = TRUE)]], " x ", colnames(labelmatrix)[col(labelmatrix)[lower.tri(labelmatrix, diag = TRUE)]])
    md.write <- matrix(0, nrow = dim(object[["D"]])[3], ncol = length(dlabels) + 1)
    colnames(md.write) <- c("Iteration", dlabels)
    md.write[, "Iteration"] <- as.numeric(dimnames(object[["D"]])[[3]])
    
    for (i in 1:dim(object[["D"]])[3]){ 
      md.write[i, -1] <- object[["D"]][, , i][lower.tri(object[["D"]][, , i], diag = TRUE)]
    }
    write.table(signif(md.write, gSIGDIG), paste0(modelname, "_D.csv"), sep = ",", row.names = FALSE)
  }
  if (!is.null(object[["Draws"]]) & writeDraws) {
    cat("Creating individual draw files, this may take a few minutes.\n")
    for (i in 1:length(object[["Draws"]])) {
      fn <- paste0("Draws_", names(object[["Draws"]])[[i]], ".csv")
      write.table(signif(object[["Draws"]][[i]], gSIGDIG), fn, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  }
  
  setwd(orig.path)
}
