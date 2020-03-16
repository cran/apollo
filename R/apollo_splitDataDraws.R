#' Splits data and draws for loading in cluster
#' 
#' Splits apollo_inputs into pieces and writes them to disk (temporary folder).
#' 
#' Internal use only. This function is called by \link{apollo_makeCluster}.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If TRUE, no information is printed to console or default output.
#' @return Character vector of file names of the pieces of apollo_inputs.
apollo_splitDataDraws <- function(apollo_inputs, silent=FALSE){
  apollo_control = apollo_inputs[["apollo_control"]]
  database       = apollo_inputs[["database"]]
  
  ### Extract useful elements
  nObs    <- nrow(database)
  indivID <- database[,apollo_control$indivID]
  namesID <- unique(indivID)
  nIndiv  <- length(namesID)
  nObsID  <- as.vector(table(indivID))
  mixing  <- apollo_control$mixing
  nCores  <- apollo_control$nCores
  
  ### Figure out how to split the database per individual
  database <- database[order(indivID),]
  indivID  <- database[,apollo_control$indivID] # update order
  namesID  <- unique(indivID)                   # update order
  apollo_inputs$database <- database            # update order
  rm(database); gc()
  if(!silent) cat('\nAttempting to split data into', nCores, 'pieces.',sep=' ')
  obj          <- ceiling(nObs/nCores)
  counter      <- 0
  currentCore  <- 1
  assignedCore <- rep(0,nObs) 
  i <- 1
  for(n in 1:nIndiv){
    assignedCore[i:(i+nObsID[n]-1)] <- currentCore
    i <- i+nObsID[n]
    counter <- counter + nObsID[n]
    if(counter>=obj & currentCore<nCores){
      currentCore <- currentCore + 1
      counter <- 0
    }
  }
  nCores <- max(assignedCore)
  if(!silent && nCores!=apollo_inputs$apollo_control$nCores) cat("\n ",nCores, ' workers (threads) will be used for estimation.', sep='')
  apollo_inputs$apollo_control$nCores <- nCores
  coreLoad <- as.vector(table(assignedCore)) 
  names(coreLoad) <- paste('worker',1:nCores,sep='_')
  if(!silent) cat('\n Number of observations per worker (thread):\n  ', 
                  paste0(coreLoad, collapse=", "), sep="")
  rm(obj, counter, currentCore, i, n)
  
  
  ### Create names of temp files
  inputPieceFile <- rep("", nCores)
  for(i in 1:nCores) inputPieceFile[i] <- tempfile()
  
  ### Identify elements to split
  # For lists, it looks at its first element
  byObs <- c()
  byInd <- c()
  asIs  <- c()
  for(e in ls(apollo_inputs)){
    E <- apollo_inputs[[e]]
    if(is.list(E)) E <- E[[1]]
    if(is.array(E)) nRows <- dim(E)[1] else nRows <- length(E)
    if(nRows==nObs) byObs <- c(byObs, e)
    if(nRows==nIndiv) byInd <- c(byInd, e)
    if(nRows!=nObs & nRows!=nIndiv) asIs <- c(asIs, e)
  }; rm(e, E, nRows)
  
  ### Create and write to disk each fragment of apollo_inputs
  if(!silent) cat("\n Writing pieces to disk")
  for(i in 1:nCores){
    L <- list()
    # Elements to be kept unchanged
    L <- apollo_inputs[ asIs ]
    # Elements to be splitted by observations
    if(length(byObs)>0) for(e in byObs){
      E <- apollo_inputs[[e]]
      r <- assignedCore==i
      if(is.list(E) & !is.data.frame(E)){
        L[[e]] <- lapply(E, apollo_keepRows, r=r)
      } else L[[e]] <- apollo_keepRows(E, r)
    }
    # Elements to be splitted by individual
    if(length(byInd)>0) for(e in byInd){
      E <- apollo_inputs[[e]]
      r <- namesID %in% indivID[assignedCore==i]
      if(is.list(E) & !is.data.frame(E)){
        L[[e]] <- lapply(E, apollo_keepRows, r=r)
      } else L[[e]] <- apollo_keepRows(E, r)
    }
    # Write to disk
    wroteOK <- tryCatch({
      saveRDS(L, file=inputPieceFile[i])
      TRUE
    }, warning=function(w) FALSE, error=function(e) FALSE)
    if(!wroteOK) stop("Apollo could not write data pieces to disk. This is necessary for multi-core processing to work.")
    if(!silent) cat(".")
  }
  
  if(!silent) cat(" Done. ", sum(gc()[,2]), 'MB of RAM in use.', sep="")
  return(inputPieceFile)
}