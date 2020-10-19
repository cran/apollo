#' Calculate DFT probabilities
#' 
#' Calculate probabilities of a Decision Field Theory (DFT) with external thresholds.
#' 
#' @param dft_settings List of settings for the DFT model. It should contain the following elements.
#'                      \itemize{
#'                        \item \strong{alternatives}: Named numeric vector. Names of alternatives and 
#'                                                     their corresponding value in \code{choiceVar}.
#'                        \item \strong{avail}: Named list of numeric vectors or scalars. Availabilities 
#'                                              of alternatives, one element per alternative. Names of 
#'                                              elements must match those in \code{alternatives}. Values 
#'                                              can be 0 or 1.
#'                        \item \strong{choiceVar}: Numeric vector. Contains choices for all observations. 
#'                                                  It will usually be a column from the database. Values 
#'                                                  are defined in \code{alternatives}.
#'                        \item \strong{attrValues}: A named list with as many elements as alternatives. 
#'                                                   Each element is itself a named list of vectors of the 
#'                                                   alternative attributes for each observation (usually a 
#'                                                   column from the database). All alternatives must have 
#'                                                   the same attributes (can be set to zero if not relevant).
#'                        \item \strong{altStart}: A named list with as many elements as alternatives. 
#'                                                 Each elment can be a scalar or vector containing the 
#'                                                 starting preference value for the alternative.  
#'                        \item \strong{attrWeights}: A named list with as many elements as attributes, 
#'                                                    or fewer. Each element is the weight of the attribute, 
#'                                                    and can be a scalar, a vector with as many elements as 
#'                                                    observations, or a matrix/array if random. They should 
#'                                                    add up to one for each observation and draw (if present), 
#'                                                    and will be re-scaled if they do not. \code{attrWeights} 
#'                                                    and \code{attrScalings} are incompatible, and they should 
#'                                                    not be both defined for an attribute. Default is 1 for 
#'                                                    all attributes.
#'                        \item \strong{attrScalings}: A named list with as many elements as attributes, 
#'                                                     or fewer. Each element is a factor that scale the 
#'                                                     attribute, and can be a scalar, a vector or a 
#'                                                     matrix/array. They do not need to add up to one 
#'                                                     for each observation. \code{attrWeights} and 
#'                                                     \code{attrScalings} are incompatible, and they 
#'                                                     should not be both defined for an attribute. 
#'                                                     Default is 1 for all attributes.
#'                        \item \strong{procPars}: A list containing the four DFT 'process parameters'
#'                          \itemize{
#'                            \item \strong{error_sd}: Numeric scalar or vector. The standard deviation of the the error term in each timestep.
#'                            \item \strong{timesteps}: Numeric scalar or vector. Number of timesteps to consider. Should be an integer bigger than 0.
#'                            \item \strong{phi1}: Numeric scalar or vector. Sensitivity.
#'                            \item \strong{phi2}: Numeric scalar or vector. Process parameter.
#'                          }
#'                       \item \strong{rows}: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                       \item \strong{componentName}: Character. Name given to model component.
#'                      }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item \code{"estimate"}: Used for model estimation.
#'                        \item \code{"prediction"}: Used for model predictions.
#'                        \item \code{"validate"}: Used for validating input.
#'                        \item \code{"zero_LL"}: Used for calculating null likelihood.
#'                        \item \code{"conditionals"}: Used for calculating conditionals.
#'                        \item \code{"output"}: Used for preparing output after model estimation.
#'                        \item \code{"raw"}: Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item \strong{\code{"estimate"}}: vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item \strong{\code{"prediction"}}: List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item \strong{\code{"validate"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"zero_LL"}}: vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item \strong{\code{"conditionals"}}: Same as \code{"estimate"}
#'           \item \strong{\code{"output"}}: Same as \code{"estimate"} but also writes summary of input data to internal Apollo log.
#'           \item \strong{\code{"raw"}}: Same as \code{"prediction"}
#'         }
#' @section References:
#' Hancock, T.; Hess, S. and Choudhury, C. (2018) Decision field theory: Improvements to current methodology and comparisons with standard choice modelling techniques. Transportation Research 107B, 18 - 40.
#' Hancock, T.; Hess, S. and Choudhury, C. (Submitted) An accumulation of preference: two alternative dynamic models for understanding transport choices.
#' Roe, R.; Busemeyer, J. and Townsend, J. (2001) Multialternative decision field theory: A dynamic connectionist model of decision making. Psychological Review 108, 370
#' @export
#' @importFrom mnormt pmnorm
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @import Rcpp
#' @useDynLib apollo, .registration=TRUE
apollo_dft = function(dft_settings,functionality){
  ### Set or extract componentName
  modelType = "DFT"
  if(is.null(dft_settings[["componentName"]])){
    dft_settings[["componentName"]] = ifelse(!is.null(dft_settings[['componentName2']]),
                                             dft_settings[['componentName2']], modelType)
    test <- functionality=="validate" && dft_settings[["componentName"]]!='model' && !apollo_inputs$silent
    if(test) apollo_print(paste0('Apollo found a model component of type ', modelType,
                                 ' without a componentName. The name was set to "',
                                 dft_settings[["componentName"]],'" by default.'))
  }
  ### Check for duplicated modelComponent name
  if(functionality=="validate"){
    apollo_modelList <- tryCatch(get("apollo_modelList", envir=parent.frame(), inherits=FALSE), error=function(e) c())
    apollo_modelList <- c(apollo_modelList, dft_settings$componentName)
    if(anyDuplicated(apollo_modelList)) stop("Duplicated componentName found (", dft_settings$componentName,
                                             "). Names must be different for each component.")
    assign("apollo_modelList", apollo_modelList, envir=parent.frame())
  }
  
  # ############################### #
  #### Load or do pre-processing ####
  # ############################### #
  # Fetch apollo_inputs
  apollo_inputs = tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE),
                           error=function(e) return( list(apollo_control=list(cpp=FALSE)) ))
  
  if( !is.null(apollo_inputs[[paste0(dft_settings$componentName, "_settings")]]) && (functionality!="preprocess") ){
    # Load dft_settings from apollo_inputs
    tmp <- apollo_inputs[[paste0(dft_settings$componentName, "_settings")]]
    # If there is no V inside the loaded dft_settings, restore the one received as argument
    if(is.null(tmp$altStart   )) tmp$altStart    <- dft_settings$altStart   
    if(is.null(tmp[['attrWeights']])) tmp$attrWeights <- dft_settings$attrWeights
    if(is.null(tmp[['attrScalings']])) tmp$attrScalings <- dft_settings$attrScalings
    if(is.null(tmp$procPars   )) tmp$procPars    <- dft_settings$procPars   
    dft_settings <- tmp
    rm(tmp)
    
  } else { 
    ### Do pre-processing
    # Do pre-processing common to most models
    dft_settings <- apollo_preprocess(inputs = dft_settings, modelType, 
                                      functionality, apollo_inputs)
    
    # Determine which likelihood to use (R or C++)
    if(apollo_inputs$apollo_control$cpp) if(!apollo_inputs$silent) apollo_print(paste0('No C++ optimisation available for ',
                                                                                       modelType))
    dft_settings$probs_DFT <- function(dft_settings, all=FALSE){
      
      ### Not really sure how this will impact DFT code?..
      # Fix choiceVar if "raw" and choiceVar==NA
      dft_settings$choiceNA = FALSE
      if(length(dft_settings$choiceVar)==1 && is.na(dft_settings$choiceVar)){
        dft_settings$choiceVar = dft_settings$alternatives[1]
        dft_settings$choiceNA  = TRUE
      }
      
      ### Extract values
      ### other parameter adjustments?...
      ### square error here? DFT calculation requires variance, but might be simpler to expect that sd is provided
      erv  = dft_settings$procPars[["error_sd"]]^2
      ts   = dft_settings$procPars[["timesteps"]]
      phi1 = dft_settings$procPars[["phi1"]]
      phi2 = dft_settings$procPars[["phi2"]]
      
      ######### DFT Probability part
      
      ### if dimension = 1 (no mixing)
      ##### then need to rearrange all lists into matrices, so that mapply can be used
      if(dft_settings$Dims==1){
        attrValuesM <- array(c(unlist(dft_settings$attrValues,use.names=F)),
                             c(dft_settings$nObs, dft_settings$nAttrs*dft_settings$nAlt))
        availM <- array(c(unlist(dft_settings$avail,use.names=F)),c(dft_settings$nObs, dft_settings$nAlt))
        altStartM <- c()
        for(i in 1:dft_settings$nAlt){
          if(length(dft_settings$altStart[[i]])==dft_settings$nObs) {
            altStartM <- c(altStartM,dft_settings$altStart[[i]]) 
          } else {
            altStartM <- c(altStartM,rep(dft_settings$altStart[[i]],dft_settings$nObs))
          }
        }
        dim(altStartM) <- c(dft_settings$nObs,dft_settings$nAlt)
        
        if(sum(lengths(dft_settings$attrWeights))!=1){
          attrWeightsM <- c()
          for(i in 1:dft_settings$nAttrs){
            if(length(dft_settings$attrWeights[[i]])==dft_settings$nObs) {
              attrWeightsM <- c(attrWeightsM,dft_settings$attrWeights[[i]]) 
            } else {
              attrWeightsM <- c(attrWeightsM,rep(dft_settings$attrWeights[[i]],dft_settings$nObs))
            }
          }
          dim(attrWeightsM) <- c(dft_settings$nObs,dft_settings$nAttrs)
        } else attrWeightsM <- array(1/dft_settings$nAttrs,c(dft_settings$nObs,dft_settings$nAttrs))
        
        ### could still have alternative and attribute specific scalings
        if(sum(lengths(dft_settings$attrScalings))!=1){
          attrScalingsM <- c()
          for(i in 1:dft_settings$nAttrs) if(is.list(dft_settings$attrScalings[[i]])) for(j in 1:dft_settings$nAlt) {
            if (length(dft_settings$attrScalings[[i]][[j]])==dft_settings$nObs){
              attrScalingsM <- c(attrScalingsM,dft_settings$attrScalings[[i]][[j]])
            } else{
              attrScalingsM <- c(attrScalingsM,rep(dft_settings$attrScalings[[i]][[j]],each=dft_settings$nObs))
            }
          } else {
            if (length(dft_settings$attrScalings[[i]])==dft_settings$nObs) {
              attrScalingsM <- c(attrScalingsM,rep(dft_settings$attrScalings[[i]],dft_settings$nAlt))
            } else {
              attrScalingsM <- c(attrScalingsM,rep(dft_settings$attrScalings[[i]],dft_settings$nAlt*dft_settings$nObs))
            }
          }
          dim(attrScalingsM) <- c(dft_settings$nObs,dft_settings$nAttrs*dft_settings$nAlt)
        } else {
          if (dft_settings$attrScalings!=1) stop("If you are not using dft_settings$attrScalings for model component \"",
                                                 dft_settings$componentName,"\", please set it to 1")
          attrScalingsM <- array(1,c(dft_settings$nObs,dft_settings$nAttrs*dft_settings$nAlt))
        }
        
        ### create a value for each choice for the process parameters if just a single value is passed in
        if(length(erv)==1) erv = rep(erv,dft_settings$nObs)
        if(length(ts)==1) ts = rep(ts,dft_settings$nObs)
        if(length(phi1)==1) phi1 = rep(phi1,dft_settings$nObs)
        if(length(phi2)==1) phi2 = rep(phi2,dft_settings$nObs)
      }
      
      
      if(dft_settings$Dims==2){
        Dim2Length = apollo_inputs$apollo_draws$interNDraws
        #### dft_settings$attrValues (all will be 1D)
        attrValuesM <- array(c(rep(c(unlist(dft_settings$attrValues,use.names=F)),each=Dim2Length)),c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs*dft_settings$nAlt))
        #### avail (all will be 1D)
        availM<-array(c(rep(c(unlist(dft_settings$avail,use.names=F)),each=Dim2Length)),c(dft_settings$nObs*Dim2Length,dft_settings$nAlt))
        #### dft_settings$altStart
        altStartM<-c()
        for(i in 1:dft_settings$nAlt){
          if(is.null(dim(dft_settings$altStart[[i]]))) {
            if(length(dft_settings$altStart[[i]])==dft_settings$nObs) {
              altStartM<-c(altStartM,rep(dft_settings$altStart[[i]],each=Dim2Length)) 
            } else {
              altStartM<-c(altStartM,rep(dft_settings$altStart[[i]],dft_settings$nObs*Dim2Length))
            }
          } else {
            ## 2D
            if(dim(dft_settings$altStart[[i]])[1]==dft_settings$nObs){
              altStartM<-c(altStartM,t(dft_settings$altStart[[i]])) 
            } else {
              altStartM<-c(altStartM,rep(dft_settings$altStart[[i]],dft_settings$nObs))
            }
          }
        }
        dim(altStartM)<-c(dft_settings$nObs*Dim2Length,dft_settings$nAlt)
        
        #### dft_settings$attrWeights
        if(sum(lengths(dft_settings$attrWeights))!=1){
          attrWeightsM<-c()
          for(i in 1:dft_settings$nAttrs){
            if(is.null(dim(dft_settings$attrWeights[[i]]))) {
              if(length(dft_settings$attrWeights[[i]])==dft_settings$nObs) {
                attrWeightsM<-c(attrWeightsM,rep(dft_settings$attrWeights[[i]],each=Dim2Length)) 
              } else {
                attrWeightsM<-c(attrWeightsM,rep(dft_settings$attrWeights[[i]],dft_settings$nObs*Dim2Length))
              }
            } else {
              ## 2D
              if(dim(dft_settings$attrWeights[[i]])[1]==dft_settings$nObs){
                attrWeightsM<-c(attrWeightsM,t(dft_settings$attrWeights[[i]])) 
              } else {
                attrWeightsM<-c(attrWeightsM,rep(dft_settings$attrWeights[[i]],dft_settings$nObs))
              }
            }
          }
          dim(attrWeightsM)<-c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs)
        } else {
          attrWeightsM<-array(1/dft_settings$nAttrs,c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs))
        }
        
        #### dft_settings$attrScalings
        ## could not be provided (=1)
        ## then for each attribute:
        ## could be attribute specific -> not mixed, some mixed, all mixed
        ## could be general -> not mixed, some mixed, all mixed
        ## for each attribute:
        if(sum(lengths(dft_settings$attrScalings))!=1){
          attrScalingsM<-c()
          for(i in 1:dft_settings$nAttrs) {
            if (is.list(dft_settings$attrScalings[[i]])) {
              ### attribute-specific
              for (j in 1:dft_settings$nAlt){
                if(is.null(dim(dft_settings$attrScalings[[i]][[j]]))){
                  ### not mixed
                  if(length(dft_settings$attrScalings[[i]][[j]])==dft_settings$nObs){
                    attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]][[j]]),each=Dim2Length))
                  } else{
                    attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]][[j]]),each=dft_settings$nObs*Dim2Length))
                  }
                } else {
                  ### is mixed
                  if(length(dft_settings$attrScalings[[i]][[j]])==dft_settings$nObs*Dim2Length){
                    attrScalingsM<-c(attrScalingsM,c(t(dft_settings$attrScalings[[i]][[j]])))
                  } else{
                    attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]][[j]]),dft_settings$nObs))
                  }
                }
              }
            } else {
              ### is not attribute-specific
              if(is.null(dim(dft_settings$attrScalings[[i]]))){
                ## not mixed
                if(length(dft_settings$attrScalings[[i]])==dft_settings$nObs) {
                  attrScalingsM<-c(attrScalingsM,rep(c(rep(c(dft_settings$attrScalings[[i]]),each=Dim2Length)),dft_settings$nAlt))
                } else {
                  attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]]),dft_settings$nAlt*Dim2Length*dft_settings$nObs))
                }
              } else{
                ## is mixed
                if(length(dft_settings$attrScalings[[i]])==dft_settings$nObs*Dim2Length) {
                  attrScalingsM<-c(attrScalingsM,rep(c(t(dft_settings$attrScalings[[i]])),dft_settings$nAlt))
                } else {
                  attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]]),dft_settings$nAlt*dft_settings$nObs))
                }
              }
            }
          }
          ### build matrix
          attrScalingsM<-array(attrScalingsM,c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs*dft_settings$nAlt))
        } else {
          attrScalingsM<-array(1,c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs))
          if (dft_settings$attrScalings!=1) stop("If you are not using dft_settings$attrScalings for model component \"",
                                                 dft_settings$componentName,"\", please set it to 1")
        }
        
        ### process parameters:
        if(is.null(dim(erv))){
          if(length(erv)==dft_settings$nObs){
            erv = rep(erv,each=Dim2Length)
          } else {
            erv = rep(erv,Dim2Length*dft_settings$nObs)
          } 
        } else {
          if(dim(erv)[1]==dft_settings$nObs){
            erv = c(t(erv))
          } else {
            erv = rep(erv,dft_settings$nObs)
          }
        }
        
        if(is.null(dim(ts))){
          if(length(ts)==dft_settings$nObs){
            ts = rep(ts,each=Dim2Length)
          } else {
            ts = rep(ts,Dim2Length*dft_settings$nObs)
          } 
        } else {
          if(dim(ts)[1]==dft_settings$nObs){
            ts = c(t(ts))
          } else {
            ts = rep(ts,dft_settings$nObs)
          }
        }
        
        if(is.null(dim(phi1))){
          if(length(phi1)==dft_settings$nObs){
            phi1 = rep(phi1,each=Dim2Length)
          } else {
            phi1 = rep(phi1,Dim2Length*dft_settings$nObs)
          } 
        } else {
          if(dim(phi1)[1]==dft_settings$nObs){
            phi1 = c(t(phi1))
          } else {
            phi1 = rep(phi1,dft_settings$nObs)
          }
        }
        
        if(is.null(dim(phi2))){
          if(length(phi2)==dft_settings$nObs){
            phi2 = rep(phi2,each=Dim2Length)
          } else {
            phi2 = rep(phi2,Dim2Length*dft_settings$nObs)
          } 
        } else {
          if(dim(phi2)[1]==dft_settings$nObs){
            phi2 = c(t(phi2))
          } else {
            phi2 = rep(phi2,dft_settings$nObs)
          }
        }
        
        
        dft_settings$choiceVar = rep(dft_settings$choiceVar,each=Dim2Length) # CAN BE DONE AT THE BEGINNING
        
      }
      
      if(dft_settings$Dims==3){
        
        Dim2Length =  apollo_inputs$apollo_draws$interNDraws # CAN BE DONE AT THE BEGINNING
        Dim3Length =  apollo_inputs$apollo_draws$intraNDraws  # CAN BE DONE AT THE BEGINNING
        
        #### dft_settings$attrValues (all will be 1D) # CAN BE DONE AT THE BEGINNING
        attrValuesM<-array(c(rep(c(unlist(dft_settings$attrValues,use.names=F)),each=Dim2Length*Dim3Length)),c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAttrs*dft_settings$nAlt))
        
        #### avail (all will be 1D) # CAN BE DONE AT THE BEGINNING
        availM<-array(c(rep(c(unlist(dft_settings$avail,use.names=F)),each=Dim2Length*Dim3Length)),c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAlt))
        
        #### dft_settings$altStart
        altStartM<-c()
        for(i in 1:dft_settings$nAlt) {
          if(is.null(dim(dft_settings$altStart[[i]]))) {
            ### 1D, no mixing:
            if(length(dft_settings$altStart[[i]])==dft_settings$nObs){
              ### diff value for different obs
              altStartM<-c(altStartM,rep(c(dft_settings$altStart[[i]]),each=Dim2Length*Dim3Length))
            } else {
              altStartM<-c(altStartM,c(rep(c(dft_settings$altStart[[i]]),dft_settings$nObs*Dim2Length*Dim3Length)))
            } 
          } else {
            ### check if 2D:
            if (length(dim(dft_settings$altStart[[i]]))==2) {
              ### 2D
              if(dim(dft_settings$altStart[[i]])[1]==dft_settings$nObs){
                altStartM<-c(altStartM,rep(t(dft_settings$altStart[[i]]),each=Dim3Length)) 
              } else {
                altStartM<-c(altStartM,rep(rep(dft_settings$altStart[[i]],dft_settings$nObs),each=Dim3Length))
              }
            } else {
              ### 3D
              if(dim(dft_settings$altStart[[i]])[1]==dft_settings$nObs){
                ### may or may not have 2nd dim
                if(dim(dft_settings$altStart[[i]])[2]==Dim2Length){
                  altStartM<-c(altStartM,c(aperm(dft_settings$altStart[[i]])))
                } else {
                  altStartM<-c(altStartM,c(t(matrix(rep(dft_settings$altStart[[i]],each=Dim2Length),ncol=Dim3Length))))
                }
              } else {
                if(dim(dft_settings$altStart[[i]])[2]==Dim2Length){
                  altStartM<-c(altStartM,rep(c(aperm(dft_settings$altStart[[i]])),dft_settings$nObs))
                } else{
                  altStartM<-c(altStartM,rep(dft_settings$altStart[[i]],dft_settings$nObs*Dim2Length))
                }
              }
            }
          }
        }
        altStartM<-array(altStartM,c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAlt))
        
        
        #### dft_settings$attrWeights
        if(sum(lengths(dft_settings$attrWeights))!=1){
          attrWeightsM<-c()
          for(i in 1:dft_settings$nAttrs) {
            if(is.null(dim(dft_settings$attrWeights[[i]]))) {
              ### 1D, no mixing:
              if(length(dft_settings$attrWeights[[i]])==dft_settings$nObs){
                ### diff value for different obs
                attrWeightsM<-c(attrWeightsM,rep(c(dft_settings$attrWeights[[i]]),each=Dim2Length*Dim3Length))
              } else {
                attrWeightsM<-c(attrWeightsM,c(rep(c(dft_settings$attrWeights[[i]]),dft_settings$nObs*Dim2Length*Dim3Length)))
              } 
            } else {
              ### check if 2D:
              if (length(dim(dft_settings$attrWeights[[i]]))==2) {
                ### 2D
                if(dim(dft_settings$attrWeights[[i]])[1]==dft_settings$nObs){
                  attrWeightsM<-c(attrWeightsM,rep(t(dft_settings$attrWeights[[i]]),each=Dim3Length)) 
                } else {
                  attrWeightsM<-c(attrWeightsM,rep(rep(dft_settings$attrWeights[[i]],dft_settings$nObs),each=Dim3Length))
                }
              } else {
                ### 3D
                if(dim(dft_settings$attrWeights[[i]])[1]==dft_settings$nObs){
                  ### may or may not have 2nd dim
                  if(dim(dft_settings$attrWeights[[i]])[2]==Dim2Length){
                    attrWeightsM<-c(attrWeightsM,c(aperm(dft_settings$attrWeights[[i]])))
                  } else {
                    attrWeightsM<-c(attrWeightsM,c(t(matrix(rep(dft_settings$attrWeights[[i]],each=Dim2Length),ncol=Dim3Length))))
                  }
                } else {
                  if(dim(dft_settings$attrWeights[[i]])[2]==Dim2Length){
                    attrWeightsM<-c(attrWeightsM,rep(c(aperm(dft_settings$attrWeights[[i]])),dft_settings$nObs))
                  } else{
                    attrWeightsM<-c(attrWeightsM,rep(dft_settings$attrWeights[[i]],dft_settings$nObs*Dim2Length))
                  }
                }
              }
            }
          }
          attrWeightsM<-array(attrWeightsM,c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAttrs))
        } else {
          attrWeightsM<-array(1/dft_settings$nAlt,c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAttrs))
        }
        
        #### dft_settings$attrScalings
        ## could not be provided (=1)
        ## then for each attribute:
        ## could be attribute specific -> not mixed, some mixed, all mixed (2 or 3 dim)
        ## could be general -> not mixed, some mixed, all mixed (2 or 3 dim)
        ## 16 cases: attribute specific, choice specific, inter, intra...
        ## for each attribute:
        if(sum(lengths(dft_settings$attrScalings))!=1){
          attrScalingsM<-c()
          if(is.list(dft_settings$attrScalings)){
            ## attribute specific
            for(i in 1:dft_settings$nAttrs) {
              if(is.null(dim(dft_settings$attrScalings[[i]]))) {
                ### 1D, no mixing:
                if(length(dft_settings$attrScalings[[i]])==dft_settings$nObs){
                  ### diff value for different obs
                  attrScalingsM<-c(attrScalingsM,rep(c(dft_settings$attrScalings[[i]]),each=Dim2Length*Dim3Length))
                } else {
                  attrScalingsM<-c(attrScalingsM,c(rep(c(dft_settings$attrScalings[[i]]),dft_settings$nObs*Dim2Length*Dim3Length)))
                } 
              } else {
                ### check if 2D:
                if (length(dim(dft_settings$attrScalings[[i]]))==2) {
                  ### 2D
                  if(dim(dft_settings$attrScalings[[i]])[1]==dft_settings$nObs){
                    attrScalingsM<-c(attrScalingsM,rep(t(dft_settings$attrScalings[[i]]),each=Dim3Length)) 
                  } else {
                    attrScalingsM<-c(attrScalingsM,rep(rep(dft_settings$attrScalings[[i]],dft_settings$nObs),each=Dim3Length))
                  }
                } else {
                  ### 3D
                  if(dim(dft_settings$attrScalings[[i]])[1]==dft_settings$nObs){
                    ### may or may not have 2nd dim
                    if(dim(dft_settings$attrScalings[[i]])[2]==Dim2Length){
                      attrScalingsM<-c(attrScalingsM,c(aperm(dft_settings$attrScalings[[i]])))
                    } else {
                      attrScalingsM<-c(attrScalingsM,c(t(matrix(rep(dft_settings$attrScalings[[i]],each=Dim2Length),ncol=Dim3Length))))
                    }
                  } else {
                    if(dim(dft_settings$attrScalings[[i]])[2]==Dim2Length){
                      attrScalingsM<-c(attrScalingsM,rep(c(aperm(dft_settings$attrScalings[[i]])),dft_settings$nObs))
                    } else{
                      attrScalingsM<-c(attrScalingsM,rep(dft_settings$attrScalings[[i]],dft_settings$nObs*Dim2Length))
                    }
                  }
                }
              }
            }
          } else {
            ## not attribute specific
            if(is.null(dim(dft_settings$attrScalings))) {
              ### 1D, no mixing:
              if(length(dft_settings$attrScalings)==dft_settings$nObs){
                ### diff value for different obs
                attrScalingsM<-c(attrScalingsM,rep(rep(c(dft_settings$attrScalings),each=Dim2Length*Dim3Length),dft_settings$nAlt))
              } else {
                attrScalingsM<-c(attrScalingsM,rep(c(rep(c(dft_settings$attrScalings),dft_settings$nObs*Dim2Length*Dim3Length)),dft_settings$nAlt))
              } 
            } else {
              ### check if 2D:
              if (length(dim(dft_settings$attrScalings))==2) {
                ### 2D
                if(dim(dft_settings$attrScalings)[1]==dft_settings$nObs){
                  attrScalingsM<-c(attrScalingsM,rep(rep(t(dft_settings$attrScalings),each=Dim3Length),dft_settings$nAlt)) 
                } else {
                  attrScalingsM<-c(attrScalingsM,rep(rep(rep(dft_settings$attrScalings,dft_settings$nObs),each=Dim3Length),dft_settings$nAlt))
                }
              } else {
                ### 3D
                if(dim(dft_settings$attrScalings)[1]==dft_settings$nObs){
                  ### may or may not have 2nd dim
                  if(dim(dft_settings$attrScalings)[2]==Dim2Length){
                    attrScalingsM<-c(attrScalingsM,rep(c(aperm(dft_settings$attrScalings)),dft_settings$nAlt))
                  } else {
                    attrScalingsM<-c(attrScalingsM,c(rep(t(matrix(rep(dft_settings$attrScalings,each=Dim2Length),ncol=Dim3Length))),dft_settings$nAlt))
                  }
                } else {
                  if(dim(dft_settings$attrScalings)[2]==Dim2Length){
                    attrScalingsM<-c(attrScalingsM,rep(rep(c(aperm(dft_settings$attrScalings)),dft_settings$nObs),dft_settings$nAlt))
                  } else{
                    attrScalingsM<-c(attrScalingsM,rep(rep(dft_settings$attrScalings,dft_settings$nObs*Dim2Length),dft_settings$nAlt))
                  }
                }
              }
            }
          }
          ### build matrix
          attrScalingsM<-array(attrScalingsM,c(dft_settings$nObs*Dim2Length*Dim3Length,dft_settings$nAttrs*dft_settings$nAlt))
        } else {
          attrScalingsM<-array(1,c(dft_settings$nObs*Dim2Length,dft_settings$nAttrs))
          if (dft_settings$attrScalings!=1) stop("If you are not using dft_settings$attrScalings for model component \"",
                                                 dft_settings$componentName,"\", please set it to 1")
        }
        
        ### process parameters:
        ### 8 cases:
        ### 3rd, 2nd, 1st dims different
        if(is.null(dim(erv))){
          if(length(erv)==dft_settings$nObs) {
            erv = rep(erv,each=Dim2Length*Dim3Length)
          } else {
            erv = rep(erv,Dim2Length*Dim3Length*dft_settings$nObs)
          }
        } else {
          if(length(dim(erv))==2){
            if(dim(erv)[1]==dft_settings$nObs){
              erv = rep(t(erv),each=Dim3Length)
            } else {
              erv = rep(rep(erv,each=Dim3Length),dft_settings$nObs) #doesn't matter here which rep first...
            }
          } else {
            ## dim ==3, 4 more cases
            if(dim(erv)[2]==Dim2Length){
              if(dim(erv)[1]==dft_settings$nObs){
                erv<-c(aperm(erv))
              } else {
                erv<-rep(c(aperm(erv)),dft_settings$nObs)
              }
            } else {
              if(dim(erv)[1]==dft_settings$nObs){
                erv<-c(t(matrix(rep(erv,each=Dim2Length),ncol=Dim3Length)))
              } else {
                erv<-rep(erv,dft_settings$nObs*Dim2Length)
              }
            }
          } 
        }
        
        if(is.null(dim(ts))){
          if(length(ts)==dft_settings$nObs) {
            ts = rep(ts,each=Dim2Length*Dim3Length)
          } else {
            ts = rep(ts,Dim2Length*Dim3Length*dft_settings$nObs)
          }
        } else {
          if(length(dim(ts))==2){
            if(dim(ts)[1]==dft_settings$nObs){
              ts = rep(t(ts),each=Dim3Length)
            } else {
              ts = rep(rep(ts,each=Dim3Length),dft_settings$nObs) #doesn't matter here which rep first...
            }
          } else {
            ## dim ==3, 4 more cases
            if(dim(ts)[2]==Dim2Length){
              if(dim(ts)[1]==dft_settings$nObs){
                ts<-c(aperm(ts))
              } else {
                ts<-rep(c(aperm(ts)),dft_settings$nObs)
              }
            } else {
              if(dim(ts)[1]==dft_settings$nObs){
                ts<-c(t(matrix(rep(ts,each=Dim2Length),ncol=Dim3Length)))
              } else {
                ts<-rep(ts,dft_settings$nObs*Dim2Length)
              }
            }
          } 
        }
        
        if(is.null(dim(phi1))){
          if(length(phi1)==dft_settings$nObs) {
            phi1 = rep(phi1,each=Dim2Length*Dim3Length)
          } else {
            phi1 = rep(phi1,Dim2Length*Dim3Length*dft_settings$nObs)
          }
        } else {
          if(length(dim(phi1))==2){
            if(dim(phi1)[1]==dft_settings$nObs){
              phi1 = rep(t(phi1),each=Dim3Length)
            } else {
              phi1 = rep(rep(phi1,each=Dim3Length),dft_settings$nObs) #doesn't matter here which rep first...
            }
          } else {
            ## dim ==3, 4 more cases
            if(dim(phi1)[2]==Dim2Length){
              if(dim(phi1)[1]==dft_settings$nObs){
                phi1<-c(aperm(phi1))
              } else {
                phi1<-rep(c(aperm(phi1)),dft_settings$nObs)
              }
            } else {
              if(dim(phi1)[1]==dft_settings$nObs){
                phi1<-c(t(matrix(rep(phi1,each=Dim2Length),ncol=Dim3Length)))
              } else {
                phi1<-rep(phi1,dft_settings$nObs*Dim2Length)
              }
            }
          } 
        }
        
        if(is.null(dim(phi2))){
          if(length(phi2)==dft_settings$nObs) {
            phi2 = rep(phi2,each=Dim2Length*Dim3Length)
          } else {
            phi2 = rep(phi2,Dim2Length*Dim3Length*dft_settings$nObs)
          }
        } else {
          if(length(dim(phi2))==2){
            if(dim(phi2)[1]==dft_settings$nObs){
              phi2 = rep(t(phi2),each=Dim3Length)
            } else {
              phi2 = rep(rep(phi2,each=Dim3Length),dft_settings$nObs) #doesn't matter here which rep first...
            }
          } else {
            ## dim ==3, 4 more cases
            if(dim(phi2)[2]==Dim2Length){
              if(dim(phi2)[1]==dft_settings$nObs){
                phi2<-c(aperm(phi2))
              } else {
                phi2<-rep(c(aperm(phi2)),dft_settings$nObs)
              }
            } else {
              if(dim(phi2)[1]==dft_settings$nObs){
                phi2<-c(t(matrix(rep(phi2,each=Dim2Length),ncol=Dim3Length)))
              } else {
                phi2<-rep(phi2,dft_settings$nObs*Dim2Length)
              }
            }
          } 
        }
        dft_settings$choiceVar = rep(dft_settings$choiceVar,each=Dim2Length*Dim3Length)
      }
      
      if(all){
        # Get probabilities for all alternatives
        P = list()
        if(dft_settings$Dims==1) tmp <- 1:dft_settings$nObs
        if(dft_settings$Dims==2) tmp <- 1:dft_settings$nObs*Dim2Length
        if(dft_settings$Dims==3) tmp <- 1:dft_settings$nObs*Dim2Length*Dim3Length
        for (j in 1:dft_settings$nAlt){
          tmp2 <- sapply(tmp, function(i) calculateDFTProbs(j, 
                                                            c(attrValuesM[i,]), 
                                                            c(availM[i,]), 
                                                            c(altStartM[i,]), 
                                                            c(attrWeightsM[i,]), 
                                                            c(attrScalingsM[i,]), 
                                                            erv[i], ts[i], phi1[i], phi2[i], dft_settings$nAlt, dft_settings$nAttrs))
          if(dft_settings$Dims==1) P[[dft_settings$altnames[j]]] = tmp2
          if(dft_settings$Dims==2) P[[dft_settings$altnames[j]]] = matrix(c(tmp2), dft_settings$nObs, Dim2Length, byrow=TRUE)
          if(dft_settings$Dims==3) P[[dft_settings$altnames[j]]] = aperm(array(c(tmp2), c(Dim3Length, Dim2Length, dft_settings$nObs)))
        }
        if(!dft_settings$choiceNA) P[["chosen"]] <- Reduce('+', mapply('*', dft_settings$Y, P, SIMPLIFY=FALSE))
      } else {
        # Get probability for chosen alternatives only
        if(dft_settings$Dims==1) tmp <- 1:dft_settings$nObs
        if(dft_settings$Dims==2) tmp <- 1:dft_settings$nObs*Dim2Length
        if(dft_settings$Dims==3) tmp <- 1:dft_settings$nObs*Dim2Length*Dim3Length
        tmp <- sapply(tmp, function(i) calculateDFTProbs(c(dft_settings$choiceVar[i]), 
                                                         c(attrValuesM[i,]), 
                                                         c(availM[i,]), 
                                                         c(altStartM[i,]), 
                                                         c(attrWeightsM[i,]), 
                                                         c(attrScalingsM[i,]),
                                                         erv[i], ts[i], phi1[i], phi2[i], dft_settings$nAlt, dft_settings$nAttrs))
        if(dft_settings$Dims==1) P = tmp
        if(dft_settings$Dims==2) P = matrix(c(tmp), dft_settings$nObs, Dim2Length, byrow=TRUE)
        if(dft_settings$Dims==3) P = aperm(array(c(tmp),c(Dim3Length,Dim2Length,dft_settings$nObs)))
      }
      return(P)
    }
    
    # Construct necessary input for gradient (including gradient of utilities)
    dft_settings$gradient <- FALSE
    if(dft_settings$gradient && !apollo_inputs$silent)  apollo_print(paste0('No analytical gradient available for', 
                                                                            modelType))
    
    
    # Return dft_settings if pre-processing
    if(functionality=="preprocess"){
      # Remove things that change from one iteration to the next
      dft_settings$altStart     <- NULL
      dft_settings$attrWeights  <- NULL
      dft_settings$attrScalings <- NULL
      dft_settings$procPars     <- NULL
      return(dft_settings)
    }
  }
  
  # ############################################ #
  #### Transform V into numeric and drop rows ####
  # ############################################ #
  
  ### Execute altStart, attrWeights, attrScalings & procPars
  test <- any(sapply(dft_settings$altStart, is.function))
  if(test) dft_settings$altStart = lapply(dft_settings$altStart, function(f) if(is.function(f)) f() else f )
  test <- any(sapply(dft_settings$attrWeights, is.function))
  if(test) dft_settings$attrWeights = lapply(dft_settings$attrWeights, function(f) if(is.function(f)) f() else f )
  test <- any(sapply(dft_settings$attrScalings, is.function))
  if(test) dft_settings$attrScalings = lapply(dft_settings$attrScalings, function(f) if(is.function(f)) f() else f )
  test <- any(sapply(dft_settings$procPars, is.function))
  if(test) dft_settings$procPars = lapply(dft_settings$procPars, function(f) if(is.function(f)) f() else f )
  
  ### Add zeros to altStarts for attributes not supplied.
  for(i in 1:dft_settings$nAlt) if(is.null(dft_settings$altStart[[dft_settings$altnames[i]]])){
    dft_settings$altStart[[dft_settings$altnames[i]]] = 0
  }
  
  ### Reorder altStart, attrWeights
  if(any(dft_settings$altnames != names(dft_settings$attrValues))) dft_settings$attrValues <- dft_settings$attrValues[dft_settings$altnames]
  if(is.list(dft_settings$attrScalings)) for(i in 1:dft_settings$nAttrs){
    test <- is.list(dft_settings$attrScalings[[i]]) && any(dft_settings$altnames != names(dft_settings$attrScalings[[i]]))
    if(test) dft_settings$attrScalings[[i]] <- dft_settings$attrScalings[[i]][dft_settings$altnames]
  }
  if(any(dft_settings$altnames != names(dft_settings$altStart))) dft_settings$altStart <- dft_settings$altStart[dft_settings$altnames]
  for (i in 1:dft_settings$nAlt){
    test <- is.list(dft_settings$attrValues[[i]]) && any(dft_settings$attrnames != names(dft_settings$attrValues[[i]]))
    if(test) dft_settings$attrValues[[i]] <- dft_settings$attrValues[[i]][dft_settings$attrnames]
  }
  
  ### Remove excluded rows from altStart, attrWeights, attrScalings & procPars
  if(any(!dft_settings$rows)){
    dft_settings$altStart    <- lapply(dft_settings$altStart   , apollo_keepRows, r=dft_settings$rows)
    dft_settings$attrWeights <- lapply(dft_settings$attrWeights, apollo_keepRows, r=dft_settings$rows)
    for(i in 1:length(dft_settings$attrScalings)){
      test <- is.numeric(dft_settings$attrScalings[[i]])
      if(test) dft_settings$attrScalings[[i]] <- apollo_keepRows(dft_settings$attrScalings[[i]], dft_settings$rows)
      test <- is.list(dft_settings$attrScalings[[i]])
      if(test) dft_settings$attrScalings[[i]] <- lapply(dft_settings$attrScalings[[i]], apollo_keepRows, r=dft_settings$rows)
    }
    dft_settings$procPars    <- lapply(dft_settings$procPars   , apollo_keepRows, r=dft_settings$rows)
  }
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if (functionality=="validate"){
    if(!apollo_inputs$apollo_control$noValidation) apollo_validate(dft_settings, modelType, functionality, apollo_inputs)
    
    if(!apollo_inputs$apollo_control$noDiagnostics) apollo_diagnostics(dft_settings, modelType, apollo_inputs)
    
    testL = dft_settings$probs_DFT(dft_settings, all=FALSE)
    if(any(!dft_settings$rows)) testL <- apollo_insertRows(testL, dft_settings$rows, 1) # insert excluded rows with value 1
    if(all(testL==0)) stop('All observations have zero probability at starting value for model component "', dft_settings$componentName,'"')
    if(any(testL==0) && !apollo_inputs$silent && apollo_inputs$apollo_control$debug) apollo_print(paste0('Some observations have zero probability at starting value for model component "', 
                                                                   dft_settings$componentName, '"'))
    return(invisible(testL))
  }
  
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # turn scalar availabilities into vectors
    for(i in 1:length(dft_settings$avail)){
      if(length(dft_settings$avail[[i]])==1) dft_settings$avail[[i]] <- rep(dft_settings$avail[[i]], dft_settings$nObs)
    }
    # number of available alts in each observation
    nAvAlt <- rowSums(matrix(unlist(dft_settings$avail), ncol = length(dft_settings$avail)))
    P = 1/nAvAlt # likelihood at zero
    if(any(!dft_settings$rows)) P <- apollo_insertRows(P, dft_settings$rows, 1)
    return(P)
  }
  
  
  # ############################################################ #
  #### functionality="estimate/prediction/conditionals/raw" ####
  # ############################################################ #
  
  if(functionality %in% c("estimate", "conditionals", "output", "components")){
    P <- dft_settings$probs_DFT(dft_settings, all=FALSE)
    if(any(!dft_settings$rows)) P <- apollo_insertRows(P, dft_settings$rows, 1) # insert excluded rows with value 1
    return(P)
  }
  
  if(functionality %in% c("prediction", "raw")){
    P <- dft_settings$probs_DFT(dft_settings, all=TRUE)
    # insert excluded rows with value 1
    if(any(!dft_settings$rows)) P <- lapply(P, apollo_insertRows, r=dft_settings$rows, val=1)
    return(P)
  }
  
  # ############ #
  #### Report ####
  # ############ #
  if(functionality=='report'){
    P <- list()
    apollo_inputs$silent <- FALSE
    P$data  <- capture.output(apollo_diagnostics(dft_settings, modelType, apollo_inputs, param=FALSE))
    P$param <- capture.output(apollo_diagnostics(dft_settings, modelType, apollo_inputs, data =FALSE))
    return(P)
  }
  
}


calculateDFTProbs<-function(choiceVar,attribs, avail, altStart, attrWeights, attrScalings , erv, ts, phi1, phi2, nAlt, nAttrs){
  
  M<-matrix(c(attribs),nAlt,nAttrs,byrow=TRUE)*attrScalings
  
  if(avail[choiceVar]==0) {P=0;return(P)}
  
  ### check for values in pars that will break the function-> set to a very low probability if so:
  if(is.na(phi2))          {P=0;return(P)} 
  if(is.na(phi1))          {P=0;return(P)} 
  if(is.na(ts))            {P=0;return(P)} 
  if(is.na(erv))           {P=0;return(P)} 
  if(is.infinite(phi2))    {P=0;return(P)} 
  if(is.infinite(phi1))    {P=0;return(P)} 
  if(is.infinite(ts))      {P=0;return(P)} 
  if(is.infinite(erv))     {P=0;return(P)} 
  if(any(is.na(M)))        {P=0;return(P)} 
  if(any(is.infinite(M)))  {P=0;return(P)} 
  if(phi2>=0.999)          {P=0;return(P)} 
  ### catch scenarios where M and altStart are basically zeros...
  ### check eigenvalues instead?.
  
  
  ### catch scenarios where phi2 is very close to zero- will break the function-> set to zero
  #### check that this is appropriate?...
  if(ts<1) ts=1
  if (phi2<0.0000001) phi2=0
  #if (phi1<0.0000001) phi1=0.0000001
  
  ### remove unavailable alternatives
  kp=(avail)*c(1:nAlt)
  nAvail = sum(avail)
  M=matrix(c(M[c(kp),]),nAvail,nAttrs)
  altStart = altStart[kp]
  
  ### adjust choice as required...
  newChoice = sum((choiceVar>=kp&(kp!=0)))
  
  ### apply DFT function
  MVN=DFTprob(M,altStart,attrWeights,newChoice, erv, ts, phi1, phi2)
  
  if(any(is.nan(MVN[[1]])))  {P=0;return(P)} 
  if(any(is.nan(MVN[[2]])))  {P=0;return(P)} 
  if(any(is.infinite(MVN[[1]])))  {P=0;return(P)} 
  if(any(is.infinite(MVN[[2]])))  {P=0;return(P)} 
  
  if(all(MVN[[2]]==0)) {P=1/sum(avail);return(P)} 
  if(any(diag(MVN[[2]]<=0))) {P=0;return(P)} 
  
  P=mnormt::pmnorm(x=c(MVN[1][[1]]),varcov=MVN[2][[1]])  
  
  return(P)
}

