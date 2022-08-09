#' Predicts using an estimated model
#' 
#' Calculates \code{apollo_probabilities} with functionality="prediction".
#' 
#' Structure of predictions are simplified before returning, e.g. list of vectors are turned into a matrix.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param prediction_settings List. Contains settings for this function. User input is required for all settings except those with a default or marked as optional. 
#'                            \itemize{
#'                              \item \strong{\code{modelComponent}}: Character. Name of component of apollo_probabilities output to calculate predictions for. Default is to predict for all components.
#'                              \item \strong{\code{nRep}}: Scalar integer. Only used for models that require simulation for prediction (e.g. MDCEV). Number of draws used to calculate prediction. Default is 100.
#'                              \item \strong{\code{runs}}: Numeric. Number of runs to use for computing confidence intervals of predictions.
#'                              \item \strong{\code{silent}}: Boolean. If TRUE, this function won't print any output to screen.
#'                              \item \strong{\code{summary}}: Boolean. If TRUE, a summary of the prediction is printed to screen. TRUE by default.
#'                            }
#' @param modelComponent \strong{Deprecated}. Same as \code{modelComponent} inside \code{prediction_settings}.
#' @return A list containing predictions for component \code{modelComponent} of the model described in \code{apollo_probabilities}.
#'         The particular shape of the prediction will depend on the model component.
#' @importFrom stats quantile
#' @export
apollo_prediction <- function(model, apollo_probabilities, apollo_inputs, prediction_settings=list(), modelComponent=NA){
  
  # Extract settings
  if(is.character(prediction_settings) && is.character(modelComponent)){
    modelComponent=prediction_settings
    prediction_settings=list(runs=1)
  }
  if(is.character(prediction_settings) && is.list(modelComponent)){
    modelComponent1=prediction_settings
    prediction_settings=modelComponent
    modelComponent=modelComponent1
  }
  if(!exists("prediction_settings")) prediction_settings=list()
  if(is.null(prediction_settings$modelComponent)){
    if(exists("modelComponent")) prediction_settings$modelComponent=modelComponent else prediction_settings$modelComponent=NA
  }
  if(is.null(prediction_settings$runs)) prediction_settings$runs=1
  if(is.null(prediction_settings$silent)) prediction_settings$silent=FALSE
  silent=prediction_settings$silent
  if(!is.null(apollo_inputs$silent) && apollo_inputs$silent) silent=TRUE
  if(is.null(prediction_settings$nRep)) prediction_settings$nRep <- 100L
  if(is.null(prediction_settings$summary)) prediction_settings$summary <- TRUE
  modelComponent = prediction_settings$modelComponent
  runs           = prediction_settings$runs
  apollo_inputs$nRep <- prediction_settings$nRep # Copy nRep into apollo_inputs
  
  # Validate input
  if(!is.null(model$apollo_control$HB) && model$apollo_control$HB && runs>1) stop("The calculation of confidence intervals for \'apollo_prediction\' is not applicables for models estimated using HB!") 
  if(!is.null(model$apollo_control$mixing) && model$apollo_control$mixing) apollo_print("Your model contains continuous random parameters. apollo_prediction will perform averaging across draws for these. For predictions at the level of individual draws, please call the apollo_probabilities function using model$estimate as the parameters, and with functionality=\"raw\".")
  HB = !is.null(model$apollo_control$HB) && model$apollo_control$HB
  if(!HB && (!is.numeric(model$estimate) | !is.vector(model$estimate) | is.null(names(model$estimate))))  stop("The \'model$estimates\' object should be a named vector for models not estimated using HB!")   
  if(HB){
    nObs  <- nrow(apollo_inputs$database)
    if(any(!(lengths(model$estimate)%in%c(1,nObs)))) stop("\nFor models estimated using HB, the \'model$estimate\' object needs to have one entry per row in the database for each random parameter, and a single value for each fixed parameter. This is not the case here, likely because the prediction data is different from the estimation data. Please restructure \'model$estimate\' accordingly!")
  }
  
  ### Warn the user in case elements in apollo_inputs are different from those in the global environment
  apollo_compareInputs(apollo_inputs)
  
  #### Calculate prediction at estimated parameter values
  apollo_randCoeff  = apollo_inputs[["apollo_randCoeff"]]
  apollo_lcPars     = apollo_inputs[["apollo_lcPars"]]
  apollo_checkArguments(apollo_probabilities,apollo_randCoeff,apollo_lcPars)
  if(!silent) apollo_print("Running predictions from model using parameter estimates...")
  apollo_beta  = model$estimate
  apollo_fixed = model$apollo_fixed
  predictions  = apollo_probabilities(apollo_beta, apollo_inputs, functionality="prediction")
  if(length(predictions)==1) singleElement = TRUE else singleElement = FALSE
  
  # Try to figure out type of model
  modelComponentType <- "default"
  if(length(predictions)==1) modelComponentType <- model$modelTypeList else if(anyNA(modelComponent)){
    tmp <- which(names(predictions)==modelComponent)
    if(length(tmp)>0) modelComponentType <- model$modelTypeList[tmp]
    rm(tmp)
  }
  if(length(modelComponentType)==0) modelComponentType <- "default"
  modelComponentType <- tolower(modelComponentType)
  
  # Keep only the requested component, if a component name is given
  if(!is.na(modelComponent)){
    if(is.null(predictions[[modelComponent]])) stop('A component named "', modelComponent, '" does not exist!')
    predictions <- predictions[modelComponent]
  }
  
  # Transform predictions into data.frame and drop components without predictions
  noPred <- c()
  for(m in 1:length(predictions)){
    M <- predictions[[m]]
    if(length(M)==1 && is.na(M)){
      apollo_print(paste0("Predictions do not exist for model component ", names(predictions)[m]))
      noPred <- c(noPred, m)
    } 
    if(is.list(M)){
      predictions[[m]] <- data.frame(ID = apollo_inputs$database[,apollo_inputs$apollo_control$indivID],
                                     Observation = apollo_inputs$database$apollo_sequence,
                                     do.call(cbind, M))
    } else {
    if(is.matrix(M) | is.vector(M)) predictions[[m]] <- data.frame(ID = apollo_inputs$database[,apollo_inputs$apollo_control$indivID],
                                                                   Observation = apollo_inputs$database$apollo_sequence, M)
    }
  }; rm(M)
  if(length(noPred)>0) predictions <- predictions[-noPred]
  if(length(predictions)==0){
    apollo_print('Sorry, no predictions to return.')
    return(NULL)
  }
  
  ### change 8 August
  if(!is.null(apollo_inputs$apollo_control$weights)){
  apollo_print(c("\n","Weights have been defined in apollo_control, and these will be applied to predictions too. If you want unweighted predictions despite your model using weights in estimation, you can replace the vector of weights by a vector of 1s"))
  }
  ### end change
  
  ### Check if there are weights (if there are, w is a numeric vector of weights)
  w <- is.character(apollo_inputs$apollo_control$weights)
  w <- w && (apollo_inputs$apollo_control$weights %in% names(apollo_inputs$database))
  if(w) w <- apollo_inputs$database[,apollo_inputs$apollo_control$weights]
  
  #### Prediction at the estimated values only (no runs) ####
  
  # If there are no runs, print summary and return
  if(runs==1){
    if(!silent && prediction_settings$summary){
      for(m in 1:length(predictions)){
        # Update modelComponentType
        test <- !is.null(model$modelTypeList) && !is.null(names(model$modelTypeList)) && !is.null(names(predictions))
        test <- test && names(predictions)[m] %in% names(model$modelTypeList)
        if(test){
          tmp <- which(names(model$modelTypeList)==names(predictions)[m])
          if(length(tmp)>1) tmp <- tmp[1]
          modelComponentType <- tolower( model$modelTypeList[tmp] )
        } else modelComponentType <- 'default'
        # Extract predictions
        M <- predictions[[m]]
        M <- M[,-(1:2),drop=FALSE] # remove ID and Observation
        # Print it a format appropriate to the model type
        if(modelComponentType %in% c('mdcev','mdcnev')){ # MDCEV
          if(!singleElement) apollo_print(paste0("Aggregated predictions (continuous consumption, discrete choices, and expenditure) at model estimates for model component: ", names(predictions)[m]))
          if(singleElement) apollo_print("Aggregated predictions (continuous consumption, discrete choices, and expenditure) at model estimates")
          K  <- as.integer(ncol(M)/6)
          Kn <- colnames(M)[1:K]
          Kn <- substr(Kn, 1, nchar(Kn)-10)
          tmp <- c('cont_mean', 'cont_sd', 'disc_mean', 'disc_sd', 'expe_mean', 'expe_sd')
          W <- cbind(colSums(M[, paste0(Kn, "_cont_mean")]), 
                     sqrt(colSums(M[, paste0(Kn, "_cont_sd")]^2)), 
                     colSums(M[, paste0(Kn, "_disc_mean")]), 
                     sqrt(colSums(M[, paste0(Kn, "_disc_sd")]^2)), 
                     colSums(M[, paste0(Kn, "_expe_mean")]), 
                     sqrt(colSums(M[, paste0(Kn, "_expe_sd")]^2)))
          colnames(W) <- tmp
          if(is.numeric(w)){
            UW <- M/w
            UW <- cbind(colSums(UW[, paste0(Kn, "_cont_mean")]), 
                        sqrt(colSums(UW[, paste0(Kn, "_cont_sd")]^2)), 
                        colSums(UW[, paste0(Kn, "_disc_mean")]), 
                        sqrt(colSums(UW[, paste0(Kn, "_disc_sd")]^2)), 
                        colSums(UW[, paste0(Kn, "_expe_mean")]), 
                        sqrt(colSums(UW[, paste0(Kn, "_expe_sd")]^2)))
            colnames(UW) <- paste0(tmp, "_unweighted")
            colnames(W)  <- paste0(tmp, "_weighted")
            W <- cbind(UW, W)
          }
          rownames(W) <- Kn
          M <- W
          rm(W)
        }
        if(modelComponentType=='normd'){ # Normal density
          if(!singleElement) apollo_print(paste0("Summary of predicted demand at model estimates for model component: ", names(predictions)[m]))
          if(singleElement) apollo_print("Summary of predicted demand at model estimates")
          M    <- unlist(M)
          tmp  <- c('min', '1stQ', 'median', 'mean', '3rdQ', 'max', 'aggregate')
          tmp2 <- quantile(M, probs=c(0, .25, .5, .5, .75, 1, 1), na.rm=TRUE)
          tmp2[c(4,7)] <- c(mean(M, na.rm=TRUE), sum(M, na.rm=TRUE))
          if(is.numeric(w)){
            M <- M/w
            tmp1 <- quantile(M, probs=c(0, .25, .5, .5, .75, 1, 1), na.rm=TRUE)
            tmp1[c(4,7)] <- c(mean(M, na.rm=TRUE), sum(M, na.rm=TRUE))
            M <- matrix(c(tmp1, tmp2), nrow=2, ncol=length(tmp), byrow=TRUE, 
                        dimnames=list(c("Un-weighted", "Weighted"), tmp))
          }
          if(!is.numeric(w)) M <- matrix(tmp2, nrow=1, ncol=length(tmp), dimnames=list(names(predictions)[m], tmp))
        }
        if(!(modelComponentType %in% c('mdcev', 'mdcnev', 'normd'))){ # Discrete choice
          if(!singleElement) apollo_print(paste0("Prediction at model estimates for model component: ", names(predictions)[m]))
          if(singleElement) apollo_print("Prediction at model estimates")
          if(tolower(colnames(M)[ncol(M)])=='chosen') M <- M[,-ncol(M)]
          if(is.numeric(w)) M <- rbind(`Un-weighted aggregate`=  colSums(M/w, na.rm=TRUE), 
                                       `Un-weighted average`  = colMeans(M/w, na.rm=TRUE), 
                                       `Weighted aggregate`   =  colSums(M  , na.rm=TRUE))
          if(!is.numeric(w)) M <- rbind(Aggregate=colSums(M, na.rm=TRUE), Average=colMeans(M, na.rm=TRUE))
        }
        ##print(M, digits=4)
        print(round(M,2))
        
        apollo_print("\n")
      }
      txt <- "The output from apollo_prediction is a "
      if(length(predictions)==1) txt <- paste0(txt, "matrix containing the predictions at the estimated values.") else {
        txt <- paste0(txt, 'list, with one element per model component. For each ', 
                      'model component, the list element is given by a matrix ',
                      'containing the predictions at the estimated values.')
      }
      apollo_print(txt)
    }
    if(is.list(predictions) && length(predictions)==1) predictions=predictions[[1]]
    return(predictions)
  }
  
  
  #### Prediction with multiple runs ####
  
  # Create confidence intervals by drawing from parameters asymptotic distributions
  if(anyNA(model$robvarcov)) stop('Cannot calculate confidence intervals if parameters covariance matrix contains NA values.')
  if(is.null(apollo_inputs$apollo_control$seed)) set.seed(13 + 1) else set.seed(apollo_inputs$apollo_control$seed + 1)
  beta_draws = mvtnorm::rmvnorm(n     = runs, 
                                mean  = apollo_beta[!(names(apollo_beta)%in%apollo_fixed)], 
                                sigma = model$robvarcov)
  ans <- vector(mode="list", length=length(predictions))
  names(ans) <- names(predictions)
  for(m in 1:length(predictions)){
    nObs <- nrow(predictions[[m]])
    nCol <- ncol(predictions[[m]]) - 2
    colN <- colnames(predictions[[m]])[-(1:2)]
    ans[[m]] <- list(at_estimates = predictions[[m]], 
                     draws = array(NA, dim=c(nObs, nCol, runs), dimnames=list(NULL, colN, NULL)) )
  }
  if(!silent) apollo_print("Running predictions across draws from the asymptotic distribution for maximum likelihood estimates.")
  for(r in 1:runs){
    if(!silent) apollo_print(paste0("Predicting for set of draws ", r, "/", runs, "..."))
    br = c(beta_draws[r,], apollo_beta[apollo_fixed])[names(apollo_beta)]
    Pr = apollo_probabilities(br, apollo_inputs, functionality="prediction")
    if(!is.list(Pr)) stop('apollo_probabilities(..., functionality="prediction") did not return a list, as expected')
    for(m in names(predictions)){
      if(is.list(Pr[[m]])) Pr[[m]] <- do.call(cbind, Pr[[m]])
      ans[[m]]$draws[,,r] <- Pr[[m]]
    }
  }; if(!silent) apollo_print("\n")
  
  # Print summary
  if(!silent && prediction_settings$summary){
    for(m in 1:length(ans)){
      # Update modelComponentType
      test <- !is.null(model$modelTypeList) && !is.null(names(model$modelTypeList)) && !is.null(names(predictions))
      test <- test && names(predictions)[m] %in% names(model$modelTypeList)
      if(test){
        tmp <- which(names(model$modelTypeList)==names(predictions)[m])
        if(length(tmp)>1) tmp <- tmp[1]
        modelComponentType <- tolower( model$modelTypeList[tmp] )
      } else modelComponentType <- 'default'
      # Print summary of prediction in appropriate format to each model type
      if(modelComponentType %in% c('mdcev','mdcnev')){
        if(!singleElement) apollo_print(paste0("Aggregated predictions (continuous consumption, discrete choices, and expenditure) for model component: ", names(predictions)[m]))
        if(singleElement) apollo_print(paste0("Aggregated predictions (continuous consumption, discrete choices, and expenditure)"))
        M <- ans[[m]]$at_estimates[,-(1:2)]
        K <- as.integer(ncol(M)/6)
        Kn <- colnames(M)[1:K]
        Kn <- substr(Kn, 1, nchar(Kn)-10)
        M <- matrix(colSums(M), nrow=K, ncol=6, dimnames=list(Kn, c('Cont.mean', 'Cont.sd', 
                                                                    'Disc.mean', 'Disc.sd',
                                                                    'Expend.mean', 'Expend.sd')))
        agg <- apply(ans[[m]]$draws, MARGIN=c(2,3), sum, na.rm=TRUE)[1:K,]
        cil <- apply(agg, MARGIN=1, quantile, probs=0.025)
        ciu <- apply(agg, MARGIN=1, quantile, probs=0.975)
        M <- cbind(M[,1], cil, ciu, M[,2:ncol(M)])
        colnames(M) <- c('Cont.mean', 'Quantile_0.025', 'Quantile_0.975', 'Cont.sd', 
                         'Disc.mean', 'Disc.sd', 'Expend.mean', 'Expend.sd')
        print(M, digits=4)
      }
      if(modelComponentType=='normd'){
        if(!singleElement) apollo_print(paste0("Summary of aggregated predicted demand for model component: ", names(predictions)[m]))
        if(singleElement) apollo_print(paste0("Summary of aggregated predicted demand"))
        # mean, median, aggregate \ at estimates, std.dev., quantile 0.025, quantile 0.975
        M <- matrix(NA, nrow=3, ncol=4, dimnames=list(c('Mean', 'Median', 'Aggregate'), 
                                                      c("At estimates", "Std.dev.", "Quantile 0.025", "Quantile 0.975")))
        tmp <- ans[[m]]$at_estimates[,-(1:2)]
        if(is.data.frame(tmp)) tmp <- unlist(M)
        M[,'At estimates'] <- c(mean(tmp, na.rm=TRUE), quantile(tmp, probs=0.5, na.rm=TRUE), sum(tmp, na.rm=TRUE))
        M[,'Std.dev.'] <- c( sd( apply(ans[[m]]$draws, MARGIN=c(2, 3), mean, na.rm=TRUE) ),
                             sd( apply(ans[[m]]$draws, MARGIN=c(2, 3), quantile, probs=0.5, na.rm=TRUE) ), 
                             sd( apply(ans[[m]]$draws, MARGIN=c(2, 3), sum, na.rm=TRUE) ) )
        M[,'Quantile 0.025'] <- c( quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), mean, na.rm=TRUE), probs=0.025),
                                   quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), quantile, probs=0.5, na.rm=TRUE), probs=0.025), 
                                   quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), sum, na.rm=TRUE), probs=0.025) )
        M[,'Quantile 0.975'] <- c( quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), mean, na.rm=TRUE), probs=0.975),
                                   quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), quantile, probs=0.5, na.rm=TRUE), probs=0.975), 
                                   quantile( apply(ans[[m]]$draws, MARGIN=c(2, 3), sum, na.rm=TRUE), probs=0.975) )
        print(M, digits=4)
      }
      if(!(modelComponentType %in% c('mdcev', 'mdcnev', 'normd'))){
        if(!singleElement) apollo_print(paste0("Aggregated prediction for model component: ", names(predictions)[m]))
        if(singleElement) apollo_print(paste0("Aggregated prediction"))
        est <- colSums(ans[[m]]$at_estimates[,-(1:2)], na.rm=TRUE)
        agg <- apply(ans[[m]]$draws, MARGIN=c(2,3), sum, na.rm=TRUE)
        std <- apply(agg, MARGIN=1, sd, na.rm=TRUE)
        cil <- apply(agg, MARGIN=1, quantile, probs=0.025, na.rm=TRUE)
        ciu <- apply(agg, MARGIN=1, quantile, probs=0.975, na.rm=TRUE)
        tmp <- cbind(est, rowMeans(agg), std, cil, ciu)
        colnames(tmp) <- c("at MLE", 'Sampled mean', "Sampled std.dev.", "Quantile 0.025", "Quantile 0.975")
        rownames(tmp) <- colnames(ans[[m]]$at_estimates)[-(1:2)]
        if(tolower(rownames(tmp)[nrow(tmp)])=='chosen') tmp <- tmp[-nrow(tmp),]
        print(tmp, digits=4)
        
        cat('\n')
        if(!singleElement) apollo_print(paste0("Average prediction for model component: ", names(predictions)[m]))
        if(singleElement) apollo_print(paste0("Average prediction"))
        est <- colMeans(ans[[m]]$at_estimates[,-(1:2)], na.rm=TRUE)
        avg <- apply(ans[[m]]$draws, MARGIN=c(2,3), mean, na.rm=TRUE)
        std <- apply(avg, MARGIN=1, sd, na.rm=TRUE)
        cil <- apply(avg, MARGIN=1, quantile, probs=0.025, na.rm=TRUE)
        ciu <- apply(avg, MARGIN=1, quantile, probs=0.975, na.rm=TRUE)
        tmp <- cbind(est, rowMeans(avg), std, cil, ciu)
        colnames(tmp) <- c("at MLE", 'Sampled mean', "Sampled std.dev.", "Quantile 0.025", "Quantile 0.975")
        rownames(tmp) <- colnames(ans[[m]]$at_estimates)[-(1:2)]
        if(tolower(rownames(tmp)[nrow(tmp)])=='chosen') tmp <- tmp[-nrow(tmp),]
        print(tmp, digits=4)
      }
      apollo_print("\n")
    }
    
    # Final message explaining format of output
    if(length(ans)==1){
      txt <- paste0('The output from apollo_prediction is a list with two elements: ',
                    'a data.frame containing the predictions at the estimated values, ', 
                    'and an array with predictions for different values of the parameters ',
                    'drawn from their asymptotic distribution.')
    } else {
      txt <- paste0('The output from apollo_prediction is a list, with one element per model ', 
                    'component. If the user asks for confidence intervals, then, for each ', 
                    'model component, a list with two elements is returned: ',
                    'a data.frame containing the predictions at the estimated values, and an array ', 
                    'with predictions for different values of the parameters drawn from their ',
                    'asymptotic distribution.')    }
    apollo_print(txt)
  }
  if(length(ans)==1) ans=ans[[1]]
  return(ans)
}