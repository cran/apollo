#' Splits data and draws for loading in cluster
#' 
#' Copies \code{apollo_inputs} as many times as cores, but each copy contains only part of \code{database} and \code{draws}.
#' 
#' Internal use only. This function is called by \link{apollo_makeCluster}.
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param silent Boolean. If TRUE, no information is printed to console or default output.
#' @return List. Each element is a copy of apollo_inputs, but with only a piece of \code{database} and \code{draws}.
apollo_splitDataDraws <- function(apollo_inputs, silent=FALSE){
  apollo_control = apollo_inputs[["apollo_control"]]
  database       = apollo_inputs[["database"]]
  draws          = apollo_inputs[["draws"]]
  
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
  if(!silent) cat('Attempting to split data into',nCores,'pieces.\n',sep=' ')
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
  if(!silent && nCores!=apollo_inputs$apollo_control$nCores) cat(nCores, ' workers (threads) will be used for estimation.\n', sep='')
  if(!silent) cat('Number of observations per worker (thread):\n')
  coreLoad <- as.vector(table(assignedCore)) 
  names(coreLoad) <- paste('worker',1:nCores,sep='_')
  if(!silent) print(coreLoad)
  rm(obj, counter, currentCore, i, n)
  if(!silent) cat(" ", sum(gc()[,2]),'MB of RAM in use before splitting.\n', sep="")
  
  
  ### Create list of data to be copied to each thread.
  # Each element is a copy of apollo_inputs with fewer elements in database and draws
  LL <- vector(mode="list", length=nCores)
  
  
  ### Copy unchanged elements to each apollo_inputs
  tmp <- names(apollo_inputs)
  tmp <- tmp[!(tmp %in% c("database", "draws"))]
  for(i in 1:nCores) LL[[i]] <- apollo_inputs[ tmp ]
  rm(tmp, i)
  
  
  ### Split draws and copy into corresponding apollo_inputs (LL)
  if(mixing){ 
    if(!silent) cat('Splitting draws')
    # function to split draws both for both cubes and matrices
    getDrawsPiece <- function(d,rowsID){
      if(length(dim(d))==3){
        nObsW  <- length(rowsID)
        nInter <- dim(d)[2]
        nIntra <- dim(d)[3]
        if(nInter==1 & nIntra==1) return(as.vector(d[rowsID,,])) 
        if(nInter >1 & nIntra==1) return(d[rowsID,,]) 
        if(nInter==1 & nIntra >1) return(array(d[rowsID,,], dim=c(nObsW,1,nIntra))) 
        if(nInter >1 & nIntra >1) return(d[rowsID,,]) 
      } else return(d[rowsID,])
    }
    environment(getDrawsPiece) <- new.env(parent=parent.env(environment(getDrawsPiece)))
    for(i in 1:nCores){
      LL[[i]][["draws"]] <- lapply(draws, getDrawsPiece, rowsID=which(assignedCore==i))
      if(!silent) cat(".")
    }
    rm(getDrawsPiece)
    if(!silent) cat(" Done. ", sum(gc()[,2]),'MB of RAM in use.\n', sep="")
  } else { 
    for(i in 1:nCores) LL[[i]][["draws"]] <- NA
  }
  
  
  ### Split database and copy to each apollo_inputs
  if(!silent) cat('Splitting database')
  for(i in 1:nCores){
    LL[[i]][["database"]] <- database[which(assignedCore==i),]
    if(!silent)cat('.')
  }
  if(!silent) cat(" Done. ", sum(gc()[,2]), 'MB of RAM in use.\n', sep="")
  
  return(LL)
}