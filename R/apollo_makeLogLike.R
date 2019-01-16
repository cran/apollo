#' Creates loglikelihood function.
#'
#' Creates loglikelihood function from the likelihood function apollo_probabilities provided by the user.
#'
#' Internal use only. Called by apollo_estimate before estimation.
#'
#' @param theta_fix_val Named numeric vector. Names an values of fixed parameters.
#' @param database data.frame Data to be used by model. Gets splitted and copied into the cluster.
#' @param apollo_probabilities Function. likelihood provided by the user.
#' @param apollo_control List of model options.
#' @param draws List of draws (numeric arrays). Default is NA. Get splitted and copied into the cluster.
#' @param apollo_randcoeff Function. Creates the random coeff or error terms by mixing draws and params. Provided by the user.
#' @param apollo_lcpars Function. Used for latent class models. Not implemented yet.
#' @param cl Cluster.
#' @param estimation_routine Character.
#' @param work_in_logs Boolean.
#' @return apollo_logLike function.
apollo_makeLogLike <- function(theta_fix_val, database, apollo_probabilities, apollo_control,
                            draws=NA, apollo_randcoeff=NA, apollo_lcpars=NA, cl=NA, estimation_routine, work_in_logs){

  nIter <- 0

  if(apollo_control$nCores==1){
    apollo_logLike <- function(theta_var, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){

      if(getNIter) return(nIter)

      theta <- c(theta_var, theta_fix_val)
      ### Function initialisation: do not change the following four commands
      b <- as.list(theta)
      if(apollo_control$mixing) b <- c(b, apollo_randcoeff(theta, draws))
      x <- database
      ### Calculate likelihood
      P <- apollo_probabilities(b, database, functionality="estimate")
      isVec <- is.vector(P[["model"]])
      isMat <- is.matrix(P[["model"]])
      isCub <- (is.array(P[["model"]]) && length(dim(P[["model"]]))==3)
      ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
      if(isCub) P = apollo_avgIntraDraws(P, apollo_control, functionality="estimate")
      ### Product across choices for same individual (comment out this line if not using panel data)
      if(apollo_control$panelData) P = apollo_panelProd(P,apollo_control,functionality="estimate", work_in_logs, database[,apollo_control$indivID])
      ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
      if(isMat || isCub) P = apollo_avgInterDraws(P, apollo_control, functionality="estimate", database[,apollo_control$indivID])
      ### Prepares output for function (do not change this line)
      P = apollo_prepareProb(P,apollo_control,functionality="estimate")

      if(work_in_logs & apollo_control$panelData) ll <- P  else ll <- log(P)
      if(sumLL) ll <- sum(ll)
      if(countIter & (estimation_routine=="bfgs")){
        if(exists('lastFuncParam')){
          lastFuncParam <- get("lastFuncParam", envir=.GlobalEnv)
          if(!anyNA(lastFuncParam) & all(lastFuncParam==theta_var)) nIter <<- nIter + 1
        }
      }
      #if(writeIter) apollo_writeTheta(theta_var, sum(ll), apollo_control$modelName)
      return(ll)
    }
  } else {

    rm(database, apollo_probabilities, draws)

    apollo_parProb <- function(th){
      ### Function initialisation: do not change the following four commands
      b <- as.list(th)
      if(apollo_control$mixing) b <- c(b, apollo_randcoeff(th, draws))
      x <- database
      ### Calculate likelihood
      P <- apollo_probabilities(b, database, functionality="estimate")
      isVec <- is.vector(P[["model"]])
      isMat <- is.matrix(P[["model"]])
      isCub <- (is.array(P[["model"]]) && length(dim(P[["model"]]))==3)
      ### Average across intra-individual draws (comment out this line if intra-individual mixing not used)
      if(isCub) P = apollo_avgIntraDraws(P, apollo_control, functionality="estimate")
      ### Product across choices for same individual (comment out this line if not using panel data)
      if(apollo_control$panelData) P = apollo_panelProd(P,apollo_control,functionality="estimate", work_in_logs, database[,apollo_control$indivID])
      ### Average across inter-individual draws (comment out this line if inter-individual mixing not used)
      if(isMat || isCub) P = apollo_avgInterDraws(P, apollo_control, functionality="estimate", database[,apollo_control$indivID])
      ### Prepares output for function (do not change this line)
      P = apollo_prepareProb(P,apollo_control,functionality="estimate")
      return(P)
    }

    environment(apollo_parProb) <- new.env(parent=parent.env(environment(apollo_parProb)))
    assign("work_in_logs" , work_in_logs , envir=environment(apollo_parProb))
    assign("apollo_control" , apollo_control , envir=environment(apollo_parProb))
    assign("apollo_randcoeff" , apollo_randcoeff , envir=environment(apollo_parProb))


    apollo_logLike <- function(theta_var, countIter=FALSE, writeIter=FALSE, sumLL=FALSE, getNIter=FALSE){
      if(getNIter) return(nIter)
      theta <- c(theta_var, theta_fix_val)
      P <- parallel::clusterCall(cl=cl, fun=apollo_parProb, th=theta)
      P <- unlist(P)
      if(work_in_logs & apollo_control$panelData) ll <- P  else ll <- log(P)
      if(countIter && (estimation_routine=="bfgs")){
        if(exists('lastFuncParam')){
          lastFuncParam <- get("lastFuncParam", envir=.GlobalEnv)
          if(!anyNA(lastFuncParam) & all(lastFuncParam==theta_var)) nIter <<- nIter + 1
        }
      }
      #if(writeIter) apollo_writeTheta(theta_var,sum(ll),apollo_control$modelName)
      if(sumLL) ll <- sum(ll)
      return(ll)
    }

  }

  return(apollo_logLike)
}
