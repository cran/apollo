#' Splits data and draws for loading in cluster
#' 
#' Splits the database and the draws (if given) in a list with as many elements as apollo_control$nCores.
#' 
#' Internal use only. This function is called by \link{apollo_makeCluster}.
#' @param apollo_control List. Contains options for the estimation. See \link{apollo_validatecontrol}.
#' @param database data.frame. Estimation data.
#' @param draws List. Draws as created by \link{apollo_makeDraws}.
#' @param silent Boolean. If TRUE, no information is printed to console or default output.
#' @return List. Each element is a list containing two elements:
#'         a fraction of the database and a fraction of the draws (if given).
apollo_splitDataDraws <- function(apollo_control, database, draws=NA, silent=FALSE){
  
  nObs    <- nrow(database)
  indivID <- database[,apollo_control$indivID]
  namesID <- unique(indivID)
  nIndiv  <- length(namesID)
  nObsID  <- as.vector(table(indivID))
  mixing  <- apollo_control$mixing
  nCores  <- apollo_control$nCores
  
  database <- database[order(indivID),]
  
  if(!silent) cat('Attempting to split data into',nCores,' pieces.\n',sep=' ')
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
  if(!silent) cat(nCores, ' workers (threads) will be used for estimation.\n', sep='')
  if(!silent) cat('Worker load (number of observations per thread):\n')
  coreLoad <- as.vector(table(assignedCore)) 
  names(coreLoad) <- paste('worker',1:nCores,sep='_')
  if(!silent) print(coreLoad)
  rm(obj, counter, currentCore, i, n)
  
  LL <- vector(mode="list", length=nCores)
  
  if(!silent){
    mbRAM <- sum(gc()[,2])
    cat('Splitting data. (',mbRAM,'Mb of RAM in use before splitting)', sep='')
  }
  if(apollo_control$mixing){ 
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
      rowsID <- which(assignedCore==i)
      LL[[i]] <- list(database=database[rowsID,],
                      draws=c(lapply(draws[-length(draws)], getDrawsPiece, rowsID=rowsID),
                              draws[[length(draws)]]))
    }
  } else { 
    for(i in 1:nCores){
      rowsID <- which(assignedCore==i)
      LL[[i]] <- list(database=database[rowsID,],
                      draws=NA)
    }
    if(!silent)cat('.')
  }
  if(!silent){
    mbRAM <- sum(gc()[,2])
    cat('\nData splitting completed. (',mbRAM,'Mb of RAM in use)\n\n', sep='')
  }
  
  return(LL)
}