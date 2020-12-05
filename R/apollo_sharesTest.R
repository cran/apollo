#' Compares predicted and observed shares
#' 
#' Prints tables comparing the shares predicted by the model with the shares observed in the data.
#' 
#' This is an auxiliary function to help guide the definition of utility functions in a choice model. 
#' By comparing the predicted and observed shares of alternatives for different categories of the data, 
#' it is possible to identify what additional explanatory variables could improve the fit of the model.
#' @param model Model object. Estimated model object as returned by function \link{apollo_estimate}.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Must receive three arguments:
#'                          \itemize{
#'                            \item \code{apollo_beta}: Named numeric vector. Names and values of model parameters.
#'                            \item \code{apollo_inputs}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \code{functionality}: Character. Can be either "estimate" (default), "prediction", "validate", "conditionals", "zero_LL", or "raw".
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param sharesTest_settings List of arguments. It must include the following.
#'                            \itemize{
#'                              \item \code{alternatives}: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                              \item \code{choiceVar}: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                              \item \code{subsamples}: Named list of boolean vectors. Each element of the list defines whether a given observation belongs to a given subsample (e.g. by sociodemographics).
#'                              \item \code{modelComponent}: Name of model component. Set to model by default.
#'                              \item \code{newAlts}: Optional list describing the new alternatives to be used by apollo_sharesTest. This should have as many elements as new alternatives, with each entry being a matrix of 0-1 entries, with one row per observation, and one column per alternative used in the model.
#'                              \item \code{newAltsOnly}: Boolean. If TRUE, results will only be printed for the 'new' alternatives defined in newAlts, not the original alternatives used in the model. Set to FALSE by default.
#'                            }
#' @return Nothing
#' @export
apollo_sharesTest=function(model, apollo_probabilities, apollo_inputs, sharesTest_settings){
  if(is.null(sharesTest_settings[["alternatives"]])) stop("The sharesTest_settings list needs to include an object called \"alternatives\"!")
  if(is.null(sharesTest_settings[["choiceVar"]])) stop("The sharesTest_settings list needs to include an object called \"choiceVar\"!")
  if(is.null(sharesTest_settings[["subsamples"]])) sharesTest_settings[["subsamples"]]=NA
  if(is.null(sharesTest_settings[["modelComponent"]])) sharesTest_settings$modelComponent="model"
  if(is.null(sharesTest_settings[["newAlts"]])) sharesTest_settings$newAlts=NULL  ### 31 Oct
  if(is.null(sharesTest_settings[["newAltsOnly"]])) sharesTest_settings$newAltsOnly=FALSE  ### 31 Oct
  
  alternatives = sharesTest_settings[["alternatives"]]
  choiceVar    = sharesTest_settings[["choiceVar"]] 
  subsamples   = sharesTest_settings[["subsamples"]]
  newAlts      = sharesTest_settings[["newAlts"]]  ### 31 Oct
  newAltsOnly  = sharesTest_settings[["newAltsOnly"]]  ### 31 Oct
  
  predictedShares = apollo_prediction(model, apollo_probabilities, apollo_inputs, 
                                      prediction_settings=list(modelComponent=sharesTest_settings$modelComponent,silent=TRUE))[[1]]
  
  ### 31 Oct
  if(!is.null(newAlts)){
   nAlts=ncol(predictedShares)-3
   nObs=nrow(predictedShares)
   if(any(lapply(newAlts,nrow)!=nObs)) stop("Some components in sharesTestSettings$newAlts do not have number of rows equal to that in the database!")
   if(any(lapply(newAlts,ncol)!=nAlts)) stop("Some components in sharesTestSettings$newAlts do not have number of columns equal to the number of alternatives in the model!")
  }
  
  #### predictedShares=predictedShares[,-ncol(predictedShares)] ### removed
  predictedShares=predictedShares[,!colnames(predictedShares)%in%c("ID","Observation","chosen")]#### NEW
  
  ### SH new lines 4/4. Sort out full sample here
  if(length(subsamples)==1 && is.na(subsamples)) subsamples=list("All data"=rep(1,length(choiceVar)))
  
  ### If there are any NA in predicted Shares, then it assumes these come from rows
  ### and removes them from the analysis
  if(anyNA(predictedShares)){
    rows <- !is.na(rowSums(predictedShares))
    choiceVar <- choiceVar[rows]
    subsamples <- lapply(subsamples, function(x) x[rows])
    predictedShares <- predictedShares[rows,]
    apollo_print("Warning: Predicted values contain NA. This could be due to using the rows setting in the model. These observations will be ommited from the analysis.")
  }
  
  xnames=colnames(predictedShares)
  predictedShares=split(predictedShares, rep(1:ncol(predictedShares), each = nrow(predictedShares)))
  names(predictedShares)=xnames
  
  categories = subsamples

  ### Check that values in 'categories' are either 0/1 or boolean
  isValid <- function(x){
    ux <- unique(x)
    if(all(ux %in% 0:1) || is.logical(ux)) return(TRUE)
    return(FALSE)
  }
  txt <- "Subsamples must be defined by logical (boolean) or dummy (0/1) variables."
  if(is.data.frame(categories) || is.list(categories)){
    if(!all(sapply(categories, isValid))) stop(txt)
  }
  if(is.array(categories)) if(!isValid(as.vector(categories))) stop(txt)
  if(!is.list(categories) && is.vector(categories)) if(!isValid(categories)) stop(txt)
  
  ### Calculate shares
  trueShares = list()
  for(i in 1:length(alternatives)) trueShares[[names(alternatives)[i]]] <- (choiceVar==alternatives[i])
  
  # SH changes 4/4/20
  #if(anyNA(categories)){
  #  categories=list()
  #  categories[["All data"]]=rep(1,length(trueShares[[1]]))
  #} else {
    if(!all(names(alternatives) %in% names(predictedShares))) stop("\nPredicted choice probabilities should be provided for all alternatives.")
    if(any(lapply(categories,sum)==0)) stop("\nSome subsamples are empty!)")
    if(!any(lapply(categories,sum)==length(trueShares[[1]])) & !any(lapply(categories,length)==1)) categories[["All data"]]=rep(1,length(trueShares[[1]]))
  #}
  
  cat("\nRunning share prediction tests\n")
  if(!newAltsOnly){  ## 31 Oct
    if(!is.null(newAlts)) cat("\nPart 1: alternatives as used in model")
  iterations=length(categories)
  for(j in 1:iterations){
    output=matrix(0,nrow=6,ncol=length(trueShares)+1)  
    colnames(output)=c(names(trueShares),"All")
    rownames(output)=c("Times chosen (data)","Times chosen (prediction)","SD prediction","Diff (prediction-data)","t-ratio","p-val")
    for(k in 1:length(trueShares)){
      temp1=as.data.frame(cbind(trueShares[[k]],categories[[j]]))
      temp2=as.data.frame(cbind(predictedShares[[k]],categories[[j]]))
      output[1,k]=colSums(temp1[temp1$V2==1,])[1]
      output[2,k]=colSums(temp2[temp2$V2==1,])[1]
      output[3,k]=sqrt(sum(temp2[temp2$V2==1,1]*(1-temp2[temp2$V2==1,1])))
      output[4,k]=output[2,k]-output[1,k]
      if(output[1,k]==output[2,k]){
        output[5,k]=NA
        output[6,k]=NA
      } else {
        output[5,k]=(output[2,k]-output[1,k])/output[3,k]
        output[6,k]=round(2*(1-stats::pnorm(abs(output[5,k]))),3)  
      }
    }
    output[1,(length(trueShares)+1)]=sum(output[1,(1:length(trueShares))])
    output[2,(length(trueShares)+1)]=sum(output[2,(1:length(trueShares))])
    output[3,(length(trueShares)+1)]=NA
    output[4,(length(trueShares)+1)]=sum(output[4,(1:length(trueShares))])
    output[5,(length(trueShares)+1)]=NA
    output[6,(length(trueShares)+1)]=NA
    cat("\nPrediction test for group: ",names(categories)[[j]]," (",sum(categories[[j]])," observations)",sep="")
    cat("\n")
    cat("\n")
    output=output[-3,]
    print(round(output,3))
    if(output[1,k+1]!=sum(categories[[j]])){
      cat("\nWarning: the totals in the final column are not equal to the number of observations!")
      cat("\nThis is normal if you're working with only a subset of possible alternatives in the columns.\n")
    }
  }}
  
  if(!is.null(newAlts)){
    if(!newAltsOnly) cat("\nPart 2: alternatives as defined in sharesTestSettings$newAlts")
    trueSharesNew=list()
    predictedSharesNew=list()
    for(s in 1:length(newAlts)){
      trueSharesNew[[s]]=0
      predictedSharesNew[[s]]=0
      for(j in 1:length(trueShares)){
        trueSharesNew[[s]]=trueSharesNew[[s]]+trueShares[[j]]*newAlts[[s]][,j]
        predictedSharesNew[[s]]=predictedSharesNew[[s]]+predictedShares[[j]]*newAlts[[s]][,j]
      }
    }
    names(trueSharesNew)=names(newAlts)
    names(predictedSharesNew)=names(newAlts)

    temp=Reduce("+",newAlts)
    for(j in 1:ncol(temp)){
      if(min(temp[,j])<1) cat("\nWarning: there are cases in your data where alternative",names(alternatives)[j],"does not form part of any of the new alternatives in sharesTestSettings$newAlts")
      if(max(temp[,j])>1) cat("\nWarning: there are cases in your data where alternative",names(alternatives)[j],"forms part of more than one of the new alternatives in sharesTestSettings$newAlts")
    }
    
    iterations=length(categories)
    for(j in 1:iterations){
      output=matrix(0,nrow=6,ncol=length(trueSharesNew)+1)  
      colnames(output)=c(names(trueSharesNew),"All")
      rownames(output)=c("Times chosen (data)","Times chosen (prediction)","SD prediction","Diff (prediction-data)","t-ratio","p-val")
      for(k in 1:length(trueSharesNew)){
        temp1=as.data.frame(cbind(trueSharesNew[[k]],categories[[j]]))
        temp2=as.data.frame(cbind(predictedSharesNew[[k]],categories[[j]]))
        output[1,k]=colSums(temp1[temp1$V2==1,])[1]
        output[2,k]=colSums(temp2[temp2$V2==1,])[1]
        output[3,k]=sqrt(sum(temp2[temp2$V2==1,1]*(1-temp2[temp2$V2==1,1])))
        output[4,k]=output[2,k]-output[1,k]
        if(output[1,k]==output[2,k]){
          output[5,k]=NA
          output[6,k]=NA
        } else {
          output[5,k]=(output[2,k]-output[1,k])/output[3,k]
          output[6,k]=round(2*(1-stats::pnorm(abs(output[5,k]))),3)  
        }
      }
      output[1,(length(trueSharesNew)+1)]=sum(output[1,(1:length(trueSharesNew))])
      output[2,(length(trueSharesNew)+1)]=sum(output[2,(1:length(trueSharesNew))])
      output[3,(length(trueSharesNew)+1)]=NA
      output[4,(length(trueSharesNew)+1)]=sum(output[4,(1:length(trueSharesNew))])
      output[5,(length(trueSharesNew)+1)]=NA
      output[6,(length(trueSharesNew)+1)]=NA
      cat("\nPrediction test for group: ",names(categories)[[j]]," (",sum(categories[[j]])," observations)",sep="")
      cat("\n")
      cat("\n")
      output=output[-3,]
      print(round(output,3))
      if(output[1,k+1]!=sum(categories[[j]])){
        cat("\nWarning: the totals in the final column are not equal to the number of observations!")
        cat("\nThis is normal if you're working with only a subset of possible alternatives in the columns.\n")
      }
    }
  
  }
}
