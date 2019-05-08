#' Calculate DFT probabilities
#' 
#' Calculate probabilities of a Decision Field Theory (DFT) with external thresholds.
#' 
#' @param dft_settings List of settings for the DFT model. It should contain the following elements.
#'                      \itemize{
#'                        \item alternatives: Named numeric vector. Names of alternatives and their corresponding value in \code{choiceVar}.
#'                        \item avail: Named list of numeric vectors or scalars. Availabilities of alternatives, one element per alternative. Names of elements must match those in \code{alternatives}. Values can be 0 or 1.
#'                        \item choiceVar: Numeric vector. Contains choices for all observations. It will usually be a column from the database. Values are defined in \code{alternatives}.
#'                        \item attrValues: A named list with as many elements as alternatives. Each element is itself a named list of vectors of the alternative attributes for each observation (usually a column from the database). All alternatives must have the same attributes (can be set to zero if not relevant).
#'                        \item altStart: A named list with as many elements as alternatives. Each elment can be a scalar or vector containing the starting preference value for the alternative.  
#'                        \item attrWeights: A named list with as many elements as attributes, or fewer. Each element is the weight of the attribute, and can be a scalar, a vector with as many elements as observations, or a matrix/array if random. They should add up to one for each observation and draw (if present), and will be re-scaled if they do not. \code{attrWeights} and \code{attrScalings} are incompatible, and they should not be both defined for an attribute. Default is 1 for all attributes.
#'                        \item attrScalings: A named list with as many elements as attributes, or fewer. Each element is a factor that scale the attribute, and can be a scalar, a vector or a matrix/array. They do not need to add up to one for each observation. \code{attrWeights} and \code{attrScalings} are incompatible, and they should not be both defined for an attribute. Default is 1 for all attributes.
#'                        \item procPars: A list containing the four DFT 'process parameters'
#'                          \itemize{
#'                            \item error_sd: Numeric scalar or vector. The standard deviation of the the error term in each timestep.
#'                            \item timesteps: Numeric scalar or vector. Number of timesteps to consider. Should be an integer bigger than 0.
#'                            \item phi1: Numeric scalar or vector. Sensitivity.
#'                            \item phi2: Numeric scalar or vector. Process parameter.
#'                          }
#'                       \item rows: Boolean vector. Consideration of rows in the likelihood calculation, FALSE to exclude. Length equal to the number of observations (nObs). Default is \code{"all"}, equivalent to \code{rep(TRUE, nObs)}.
#'                      }
#' @param functionality Character. Can take different values depending on desired output.
#'                      \itemize{
#'                        \item "estimate": Used for model estimation.
#'                        \item "prediction": Used for model predictions.
#'                        \item "validate": Used for validating input.
#'                        \item "zero_LL": Used for calculating null likelihood.
#'                        \item "conditionals": Used for calculating conditionals.
#'                        \item "output": Used for preparing output after model estimation.
#'                        \item "raw": Used for debugging.
#'                      }
#' @return The returned object depends on the value of argument \code{functionality} as follows.
#'         \itemize{
#'           \item "estimate": vector/matrix/array. Returns the probabilities for the chosen alternative for each observation.
#'           \item "prediction": List of vectors/matrices/arrays. Returns a list with the probabilities for all alternatives, with an extra element for the chosen alternative probability.
#'           \item "validate": Boolean. Returns TRUE if all tests are passed.
#'           \item "zero_LL": vector/matrix/array. Returns the probability of the chosen alternative when all parameters are zero.
#'           \item "conditionals": Same as "prediction".
#'           \item "output": Same as "estimate" but also writes summary of choices into temporary file (later read by \code{apollo_modelOutput}).
#'           \item "raw": Same as "prediction".
#'         }
#' @section References:
#' Hancock, T.; Hess, S. and Choudhury, C. (2018) Decision field theory: Improvements to current methodology and comparisons with standard choice modelling techniques. Transportation Research 107B, 18 - 40.
#' Hancock, T.; Hess, S. and Choudhury, C. (Submitted) An accumulation of preference: two alternative dynamic models for understanding transport choices.
#' Roe, R.; Busemeyer, J. and Townsend, J. (2001) Multialternative decision field theory: A dynamic connectionist model of decision making. Psychological Review 108, 370
#' @export
#' @importFrom mnormt pmnorm
#' @importFrom stats setNames
#' @import Rcpp
#' @useDynLib apollo, .registration=TRUE
apollo_dft = function(dft_settings,functionality){
  
  if(!exists("DFTprob")) {
    Rcpp::sourceCpp("apolloDFT.cpp")
    cat("Loading DFT c++ code, this may take a while")
  }
  
  if(is.null(dft_settings[["alternatives"]])) stop("The \"dft_settings\" list needs to include an object called \"alternatives\"!")
  if(is.null(dft_settings[["avail"]])) stop("The \"dft_settings\" list needs to include an object called \"avail\"!")
  if(is.null(dft_settings[["choiceVar"]])) stop("The \"dft_settings\" list needs to include an object called \"choiceVar\"!")
  if(is.null(dft_settings[["attrValues"]])) stop("The \"dft_settings\" list needs to include an object called \"attrValues\"!")
  if(is.null(dft_settings[["altStart"]])) stop("The \"dft_settings\" list needs to include an object called \"altStart\"!")
  if(is.null(dft_settings[["attrWeights"]])) stop("The \"dft_settings\" list needs to include an object called \"attrWeights\"!")
  if(is.null(dft_settings[["attrScalings"]])) stop("The \"dft_settings\" list needs to include an object called \"attrScalings\"!")
  if(is.null(dft_settings[["procPars"]])) stop("The \"dft_settings\" list needs to include an object called \"procPars\"!")
  if(is.null(dft_settings[["rows"]])) dft_settings[["rows"]]="all"
  
  ###### get values from DFT settings
  alternatives=dft_settings[["alternatives"]]
  avail       =dft_settings[["avail"]]
  choiceVar   =dft_settings[["choiceVar"]]
  attrValues  =dft_settings[["attrValues"]]
  altStart    =dft_settings[["altStart"]]
  attrWeights =dft_settings[["attrWeights"]]
  attrScalings=dft_settings[["attrScalings"]]
  procPars    =dft_settings[["procPars"]]
  if(is.null(dft_settings[["rows"]])) dft_settings[["rows"]]="all"
  rows        =dft_settings[["rows"]]
  
  apollo_control <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_control,
                             error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
  
  apollo_draws   <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_draws,
                             error = function(e) {
                               Dim <- function(x,d){
                                 if(is.vector(x)){ if(d==1) return(length(x)) else return(0) }
                                 if(is.array(x)){ if(d>0 && d <= length(dim(x))) return(dim(x)[d]) else return(0) }
                                 return(NA)
                               }
                               interNDraws <- max(sapply(attrWeights, Dim, d=2),
                                                  sapply(attrScalings, Dim, d=2))
                               intraNDraws <- max(sapply(attrWeights, Dim, d=3),
                                                  sapply(attrScalings, Dim, d=3))
                               tmp <- list(interNDraws=interNDraws,
                                           intraNDraws=intraNDraws)
                               return(tmp)
                             })
  
  # ############################## #
  #### functionality="validate" ####
  # ############################## #
  
  if (functionality=="validate"){
    
    #############################################################
    ## check that all four process parameter values are provided
    if(is.null(procPars[["error_sd"]])) stop("The \"procPars\" list needs to include an object called \"error_sd\"!")
    if(is.null(procPars[["timesteps"]])) stop("The \"procPars\" list needs to include an object called \"timesteps\"!")
    if(is.null(procPars[["phi1"]])) stop("The \"procPars\" list needs to include an object called \"phi1\"!")
    if(is.null(procPars[["phi2"]])) stop("The \"procPars\" list needs to include an object called \"phi2\"!")
    
    # Store useful values
    apollo_control <- tryCatch(get("apollo_inputs", parent.frame(), inherits=FALSE)$apollo_control,
                               error = function(e) list(noValidation=FALSE, noDiagnostics=FALSE))
    nObs  <- length(choiceVar)
    nAlts <- length(altStart)
    nAttrs <- max(length(attrWeights),length(attrScalings))
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    warn1 <- FALSE
    warn2 <- FALSE
    warn3 <- FALSE
    
    # check rows statement
    if(length(rows)!=length(choiceVar)) stop("The argument \"rows\" needs to either be \"all\" or a vector of length equal to the number of rows in the data!")
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), altnames)
    }
    
    # Reorder avail to match altnames order, if necessary
    if(all(altnames != names(avail))) avail <- avail[altnames]
    
    ## increase length of avail if necessary
    for (i in 1:nAlts){
      if(length(avail[[i]])==1) avail[[i]]=rep(c(avail[[i]]),nObs)
    }
    
    #### check that only one of attrWeights, attrScalings is supplied
    attrWeights=dft_settings[["attrWeights"]]
    attrScalings=dft_settings[["attrScalings"]]
    s1=sum(lengths(attrWeights))
    s2=sum(lengths(attrScalings))
    
    if(s1>1&s2>1) stop("Please set one of attrWeights or attrScalings to 1")
    if(s1>1) attrnames=names(attrWeights) else attrnames=names(attrScalings)
    
    # Check that the elements of attrValues match altnames
    if(any(!names(attrValues)%in%altnames)) stop("The \"attrValues\" attribute names do not match those given in \"alternatives\"!") 
    
    # Give warning message if any of the elements in attrValues are not in attrnames
    for (i in 1:nAlts) if(any(!names(attrValues[[i]])%in%attrnames)) warn1 <- TRUE 
    
    # Add zero for attrValues not supplied (functionality later will expand it to the correct dimension)
    for (i in 1:nAlts){
      for (j in 1:nAttrs){
        if(is.null(attrValues[[altnames[i]]][[attrnames[j]]])) attrValues[[altnames[i]]][[attrnames[j]]]=0
      }
    }
    
    # Give warning message if any of the elements in altStarts are not in altnames
    if(any(!names(altStart)%in%altnames)) warn2 <- TRUE 
    
    # Check altStart is a list:
    if(!is.list(altStart)) {
      altStart=list()
      warn3 <- TRUE
    }
    
    # Add zero for altStarts that are not supplied
    for(i in 1:nAlts){
      if(is.null(altStart[[altnames[i]]])) altStart[[altnames[i]]] = 0 
    }
    
    # Reorder attrValues, attrScalings (components) ,altStart to match altnames order, if necessary
    if(any(altnames != names(attrValues))) attrValues <- attrValues[altnames]
    
    if(is.list(attrScalings)){
      for (i in 1:nAttrs){
        if(is.list(attrScalings[[i]])) if(any(altnames != names(attrScalings[[i]]))) attrScalings[[i]] <- attrScalings[[i]][altnames]
      }
    } 
    
    if(any(altnames != names(altStart))) altStart <- altStart[altnames]
    
    # Reorder attrValues such that attribute order matches with attrScalings or attrWeights
    for (i in 1:nAlts){
      if(is.list(attrValues[[i]])) if(any(attrnames != names(attrValues[[i]]))) attrValues[[i]] <- attrValues[[i]][attrnames]
    }
    
    ### check which bits stay above, and which go into the section below
    if(apollo_control$noValidation==FALSE){
      # Check that there are at least two alternatives
      if(nAlts<2) stop("DFT requires at least two alternatives")
      
      # Check that choice vector is not empty
      if(nObs==0) stop("No choices to model")
      
      # Check that labels in choice match those in the utilities and availabilities
      choiceLabs <- unique(choiceVar)
      if(!all(altnames %in% names(avail))) stop("Alternative labels in \"altnames\" do not match those in \"avail\".")
      
      # Check that there are no values in the choice column for undefined alternatives
      if(!all(choiceLabs %in% altcodes)) stop("Value in choice column that is not included in altcodes.")
      
      # check that nothing unavailable is chosen
      chosenunavail=0
      j=1
      while(j <= length(altnames)){
        if(sum((choiceVar==altcodes[j])*(avail[[j]]==0)*rows)) chosenunavail=1
        j=j+1
      }
      if(chosenunavail==1) stop("Some alternative(s) chosen despite being listed as unavailable\n")
      
      # check that all availabilities are either 0 or 1
      for(i in 1:length(avail)) if( !all(unique(avail[[i]]) %in% 0:1) ) stop("Some availability values are not 0 or 1.")
      
      ###### further checks for dims of objects?
      
      #cat("\nAll checks passed for DFT model component\n")
      
    }
    
    if(apollo_control$noDiagnostics==FALSE){
      #if(avail_set==TRUE) warning("Availability not provided to 'apollo_dft' (or some elements are NA). Full availability assumed.")
      
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j] = sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      rownames(choicematrix) = c("Times available","Times chosen","Percentage of choice overall","Percentage of choice when available")
      colnames(choicematrix) = altnames
      #cat('Overview of choices for DFT model component:\n')
      #print(round(choicematrix,2))
      #cat("\n")
      #if(any(choicematrix[4,]==0)) cat("Warning: some alternatives are never chosen in your data!\n")
      #if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                           "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      if(warn1) content[[length(content)+1]] <- paste0("Notice: Not all of the attributes given in \"attrValues\" are \n",
                                                       paste0(rep(" ",8), collapse=""), "named in \"attrScalings\" or \"attrWeights\". These will\n", 
                                                       paste0(rep(" ",8), collapse=""), "consequently be ignored.")
      if(warn2) content[[length(content)+1]] <- paste0("Notice: Not all of the alternatives given in \"altStart\" are\n",
                                                       paste0(rep(" ",8), collapse=""),"named in \"alternatives\". These will consequently be ignored.")
      if(warn3) content[[length(content)+1]] <- paste0("Notice: A list was not supplied for \"altStart\".\n",
                                                       paste0(rep(" ",8), collapse=""),"Starting values for all alternatives will be set to zero.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for DFT model component:", content, apolloLog)
    }
    
    
    
    
    
    
    return(invisible(TRUE))
    
  }
  
  
  ####################################################
  ###### get the dimensionality of the lists
  Dims = 1
  if(is.null(apollo_control$mixing)==FALSE){
    if(apollo_control$mixing==TRUE){
      if(apollo_draws$interNDraws!=0) Dims = 2
      if(apollo_draws$intraNDraws!=0) Dims = 3
    }
  } 
  
  
  # ############################## #
  #### functionality="zero_LL" ####
  # ############################## #
  
  if(functionality=="zero_LL"){
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(alternatives)
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- vector(mode="list", length=nAlts)
      avail <- lapply(avail, function(a) 1)
      names(avail) <- altnames
    }
    
    if(!anyNA(avail)) if(all(altnames != names(avail))) avail <- avail[altnames]
    for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs) # turn scalar availabilities into vectors
    nAvAlt <- rowSums(matrix(unlist(avail), ncol = length(avail))) # number of available alts in each observation
    P = 1/nAvAlt # likelihood at zero
    P[!rows] = 1
    return(P)
  }
  
  
  
  # ############################################################ #
  #### functionality="estimate/prediction/conditionals/raw" ####
  # ############################################################ #
  
  
  if(functionality %in% c("estimate","prediction","conditionals","raw","output")){
    
    s1=sum(lengths(attrWeights))
    s2=sum(lengths(attrScalings))
    if(s1>1) attrnames=names(attrWeights) else attrnames=names(attrScalings)
    
    
    ### Not really sure how this will impact DFT code?..
    # Fix choiceVar if "raw" and choiceVar==NA
    choiceNA = FALSE
    if(length(choiceVar)==1 && is.na(choiceVar)){
      choiceVar = alternatives[1]
      choiceNA = TRUE
    }
    
    ###### get process parameter values
    ### other parameter adjustments?...
    ### square error here? DFT calculation requires variance, but might be simpler to expect that sd is provided
    erv=procPars[["error_sd"]]^2
    ts=procPars[["timesteps"]]
    phi1=procPars[["phi1"]]
    phi2=procPars[["phi2"]]
    
    # Store useful values
    nObs  <- length(choiceVar)
    nAlts <- length(altStart)
    nAttrs <- max(length(attrWeights),length(attrScalings))
    avail_set <- FALSE
    altnames=names(alternatives)
    altcodes=alternatives
    if(length(rows)==1 && rows=="all") rows=rep(TRUE,length(choiceVar))
    choiceVar[!rows]=alternatives[1]
    
    # Create availability if necessary
    if(!is.list(avail)){
      avail_set <- TRUE
      avail <- setNames(as.list(rep(1,nAlts)), altnames)
    }
    
    ### Reorder avail to match altnames order, if necessary
    if(all(altnames != names(avail))) avail <- avail[altnames]
    
    ## increase length of avail if necessary
    for (i in 1:nAlts){
      if(length(avail[[i]])==1) avail[[i]]=rep(c(avail[[i]]),nObs)
    }
    
    # Check for additional warnings
    # Give warning message if any of the elements in attrValues are not in attrnames
    for (i in 1:nAlts) if(any(!names(attrValues[[i]])%in%attrnames)) warn1 <- TRUE else warn1 <- FALSE
    # Give warning message if any of the elements in altStarts are not in altnames
    if(any(!names(altStart)%in%altnames)) warn2 <- TRUE  else warn2 <- FALSE
    # Check altStart is a list:
    if(!is.list(altStart)) warn3 <- TRUE else warn3 <- FALSE
    
    # Add zero for attrValues not supplied (functionality later will expand it to the correct dimension)
    for (i in 1:nAlts){
      for (j in 1:nAttrs){
        if(is.null(attrValues[[altnames[i]]][[attrnames[j]]])) attrValues[[altnames[i]]][[attrnames[j]]]=0
      }
    }
    
    # Add zero for altStarts that are not supplied
    for(i in 1:nAlts){
      if(is.null(altStart[[altnames[i]]])) altStart[[altnames[i]]] = 0 
    }
    
    # Reorder attrValues, attrScalings (components) ,altStart to match altnames order, if necessary
    if(any(altnames != names(attrValues))) attrValues <- attrValues[altnames]
    
    if(is.list(attrScalings)){
      for (i in 1:nAttrs){
        if(is.list(attrScalings[[i]])) if(any(altnames != names(attrScalings[[i]]))) attrScalings[[i]] <- attrScalings[[i]][altnames]
      }
    } 
    
    if(any(altnames != names(altStart))) altStart <- altStart[altnames]
    
    # Reorder attrValues such that attribute order matches with attrScalings or attrWeights
    for (i in 1:nAlts){
      if(is.list(attrValues[[i]])) if(any(attrnames != names(attrValues[[i]]))) attrValues[[i]] <- attrValues[[i]][attrnames]
    }
    
    ## check that each alternative attribute is of length nObs and update if neccessary (will extend zeros, for examples)
    for (i in 1:nAttrs){
      for (j in 1:nAlts){
        if(length(attrValues[[j]][[i]])==1) attrValues[[j]][[i]] = rep(attrValues[[j]][[i]],nObs)
      }
    }
    
    
    
    #############################################
    ######### DFT Probability part
    
    ### if dimension = 1 (no mixing)
    ##### then need to rearrange all lists into matrices, so that mapply can be used
    if (Dims==1){
      
      attrValuesM<-array(c(unlist(attrValues,use.names=F)),c(nObs,nAttrs*nAlts))
      
      availM<-array(c(unlist(avail,use.names=F)),c(nObs,nAlts))
      
      altStartM<-c()
      for(i in 1:nAlts){
        if(length(altStart[[i]])==nObs) {
          altStartM<-c(altStartM,altStart[[i]]) 
        } else {
          altStartM<-c(altStartM,rep(altStart[[i]],nObs))
        }
      }
      dim(altStartM)<-c(nObs,nAlts)
      
      if(sum(lengths(attrWeights))!=1){
        attrWeightsM<-c()
        for(i in 1:nAttrs){
          if(length(attrWeights[[i]])==nObs) {
            attrWeightsM<-c(attrWeightsM,attrWeights[[i]]) 
          } else {
            attrWeightsM<-c(attrWeightsM,rep(attrWeights[[i]],nObs))
          }
        }
        dim(attrWeightsM)<-c(nObs,nAttrs)
      } else {
        attrWeightsM<-array(1/nAttrs,c(nObs,nAttrs))
      }
      
      ### could still have alternative and attribute specific scalings
      if(sum(lengths(attrScalings))!=1){
        attrScalingsM<-c()
        for(i in 1:nAttrs) if(is.list(attrScalings[[i]])) {
          for (j in 1:nAlts){
            if (length(attrScalings[[i]][[j]])==nObs){
              attrScalingsM<-c(attrScalingsM,attrScalings[[i]][[j]])
            } else{
              attrScalingsM<-c(attrScalingsM,rep(attrScalings[[i]][[j]],each=nObs))
            }
          }
        } else {
          if (length(attrScalings[[i]])==nObs) {
            attrScalingsM<-c(attrScalingsM,rep(attrScalings[[i]],nAlts))
          } else {
            attrScalingsM<-c(attrScalingsM,rep(attrScalings[[i]],nAlts*nObs))
          }
        }
        dim(attrScalingsM)<-c(nObs,nAttrs*nAlts)
      } else {
        if (attrScalings!=1) warning("If you are not using attrScalings, please set it to 1")
        attrScalingsM<-array(1,c(nObs,nAttrs*nAlts))
      }
      
      ### create a value for each choice for the process parameters if just a single value is passed in
      if(length(erv)==1) erv = rep(erv,nObs)
      if(length(ts)==1) ts = rep(ts,nObs)
      if(length(phi1)==1) phi1 = rep(phi1,nObs)
      if(length(phi2)==1) phi2 = rep(phi2,nObs)
      
    }
    
    
    if(Dims==2){
      Dim2Length =  apollo_draws$interNDraws
      
      #### attrValues (all will be 1D)
      attrValuesM<-array(c(rep(c(unlist(attrValues,use.names=F)),each=Dim2Length)),c(nObs*Dim2Length,nAttrs*nAlts))
      
      #### avail (all will be 1D)
      availM<-array(c(rep(c(unlist(avail,use.names=F)),each=Dim2Length)),c(nObs*Dim2Length,nAlts))
      
      #### altStart
      altStartM<-c()
      for(i in 1:nAlts){
        if(is.null(dim(altStart[[i]]))) {
          if(length(altStart[[i]])==nObs) {
            altStartM<-c(altStartM,rep(altStart[[i]],each=Dim2Length)) 
          } else {
            altStartM<-c(altStartM,rep(altStart[[i]],nObs*Dim2Length))
          }
        } else {
          ## 2D
          if(dim(altStart[[i]])[1]==nObs){
            altStartM<-c(altStartM,t(altStart[[i]])) 
          } else {
            altStartM<-c(altStartM,rep(altStart[[i]],nObs))
          }
        }
      }
      dim(altStartM)<-c(nObs*Dim2Length,nAlts)
      
      #### attrWeights
      if(sum(lengths(attrWeights))!=1){
        attrWeightsM<-c()
        for(i in 1:nAttrs){
          if(is.null(dim(attrWeights[[i]]))) {
            if(length(attrWeights[[i]])==nObs) {
              attrWeightsM<-c(attrWeightsM,rep(attrWeights[[i]],each=Dim2Length)) 
            } else {
              attrWeightsM<-c(attrWeightsM,rep(attrWeights[[i]],nObs*Dim2Length))
            }
          } else {
            ## 2D
            if(dim(attrWeights[[i]])[1]==nObs){
              attrWeightsM<-c(attrWeightsM,t(attrWeights[[i]])) 
            } else {
              attrWeightsM<-c(attrWeightsM,rep(attrWeights[[i]],nObs))
            }
          }
        }
        dim(attrWeightsM)<-c(nObs*Dim2Length,nAttrs)
      } else {
        attrWeightsM<-array(1/nAttrs,c(nObs*Dim2Length,nAttrs))
      }
      
      #### attrScalings
      ## could not be provided (=1)
      ## then for each attribute:
      ## could be attribute specific -> not mixed, some mixed, all mixed
      ## could be general -> not mixed, some mixed, all mixed
      ## for each attribute:
      if(sum(lengths(attrScalings))!=1){
        attrScalingsM<-c()
        for(i in 1:nAttrs) {
          if (is.list(attrScalings[[i]])) {
            ### attribute-specific
            for (j in 1:nAlts){
              if(is.null(dim(attrScalings[[i]][[j]]))){
                ### not mixed
                if(length(attrScalings[[i]][[j]])==nObs){
                  attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]][[j]]),each=Dim2Length))
                } else{
                  attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]][[j]]),each=nObs*Dim2Length))
                }
              } else {
                ### is mixed
                if(length(attrScalings[[i]][[j]])==nObs*Dim2Length){
                  attrScalingsM<-c(attrScalingsM,c(t(attrScalings[[i]][[j]])))
                } else{
                  attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]][[j]]),nObs))
                }
              }
            }
          } else {
            ### is not attribute-specific
            if(is.null(dim(attrScalings[[i]]))){
              ## not mixed
              if(length(attrScalings[[i]])==nObs) {
                attrScalingsM<-c(attrScalingsM,rep(c(rep(c(attrScalings[[i]]),each=Dim2Length)),nAlts))
              } else {
                attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]]),nAlts*Dim2Length*nObs))
              }
            } else{
              ## is mixed
              if(length(attrScalings[[i]])==nObs*Dim2Length) {
                attrScalingsM<-c(attrScalingsM,rep(c(t(attrScalings[[i]])),nAlts))
              } else {
                attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]]),nAlts*nObs))
              }
            }
          }
        }
        ### build matrix
        attrScalingsM<-array(attrScalingsM,c(nObs*Dim2Length,nAttrs*nAlts))
      } else {
        attrScalingsM<-array(1,c(nObs*Dim2Length,nAttrs))
        if (attrScalings!=1) warning("If you are not using attrScalings, please set it to 1")
      }
      
      ### process parameters:
      if(is.null(dim(erv))){
        if(length(erv)==nObs){
          erv = rep(erv,each=Dim2Length)
        } else {
          erv = rep(erv,Dim2Length*nObs)
        } 
      } else {
        if(dim(erv)[1]==nObs){
          erv = c(t(erv))
        } else {
          erv = rep(erv,nObs)
        }
      }
      
      if(is.null(dim(ts))){
        if(length(ts)==nObs){
          ts = rep(ts,each=Dim2Length)
        } else {
          ts = rep(ts,Dim2Length*nObs)
        } 
      } else {
        if(dim(ts)[1]==nObs){
          ts = c(t(ts))
        } else {
          ts = rep(ts,nObs)
        }
      }
      
      if(is.null(dim(phi1))){
        if(length(phi1)==nObs){
          phi1 = rep(phi1,each=Dim2Length)
        } else {
          phi1 = rep(phi1,Dim2Length*nObs)
        } 
      } else {
        if(dim(phi1)[1]==nObs){
          phi1 = c(t(phi1))
        } else {
          phi1 = rep(phi1,nObs)
        }
      }
      
      if(is.null(dim(phi2))){
        if(length(phi2)==nObs){
          phi2 = rep(phi2,each=Dim2Length)
        } else {
          phi2 = rep(phi2,Dim2Length*nObs)
        } 
      } else {
        if(dim(phi2)[1]==nObs){
          phi2 = c(t(phi2))
        } else {
          phi2 = rep(phi2,nObs)
        }
      }
      
      
      choiceVar = rep(choiceVar,each=Dim2Length)
      
    }
    
    if(Dims==3){
      
      Dim2Length =  apollo_draws$interNDraws
      Dim3Length =  apollo_draws$intraNDraws 
      
      #### attrValues (all will be 1D)
      attrValuesM<-array(c(rep(c(unlist(attrValues,use.names=F)),each=Dim2Length*Dim3Length)),c(nObs*Dim2Length*Dim3Length,nAttrs*nAlts))
      
      #### avail (all will be 1D)
      availM<-array(c(rep(c(unlist(avail,use.names=F)),each=Dim2Length*Dim3Length)),c(nObs*Dim2Length*Dim3Length,nAlts))
      
      #### altStart
      altStartM<-c()
      for(i in 1:nAlts) {
        if(is.null(dim(altStart[[i]]))) {
          ### 1D, no mixing:
          if(length(altStart[[i]])==nObs){
            ### diff value for different obs
            altStartM<-c(altStartM,rep(c(altStart[[i]]),each=Dim2Length*Dim3Length))
          } else {
            altStartM<-c(altStartM,c(rep(c(altStart[[i]]),nObs*Dim2Length*Dim3Length)))
          } 
        } else {
          ### check if 2D:
          if (length(dim(altStart[[i]]))==2) {
            ### 2D
            if(dim(altStart[[i]])[1]==nObs){
              altStartM<-c(altStartM,rep(t(altStart[[i]]),each=Dim3Length)) 
            } else {
              altStartM<-c(altStartM,rep(rep(altStart[[i]],nObs),each=Dim3Length))
            }
          } else {
            ### 3D
            if(dim(altStart[[i]])[1]==nObs){
              ### may or may not have 2nd dim
              if(dim(altStart[[i]])[2]==Dim2Length){
                altStartM<-c(altStartM,c(aperm(altStart[[i]])))
              } else {
                altStartM<-c(altStartM,c(t(matrix(rep(altStart[[i]],each=Dim2Length),ncol=Dim3Length))))
              }
            } else {
              if(dim(altStart[[i]])[2]==Dim2Length){
                altStartM<-c(altStartM,rep(c(aperm(altStart[[i]])),nObs))
              } else{
                altStartM<-c(altStartM,rep(altStart[[i]],nObs*Dim2Length))
              }
            }
          }
        }
      }
      altStartM<-array(altStartM,c(nObs*Dim2Length*Dim3Length,nAlts))
      
      
      #### attrWeights
      if(sum(lengths(attrWeights))!=1){
        attrWeightsM<-c()
        for(i in 1:nAttrs) {
          if(is.null(dim(attrWeights[[i]]))) {
            ### 1D, no mixing:
            if(length(attrWeights[[i]])==nObs){
              ### diff value for different obs
              attrWeightsM<-c(attrWeightsM,rep(c(attrWeights[[i]]),each=Dim2Length*Dim3Length))
            } else {
              attrWeightsM<-c(attrWeightsM,c(rep(c(attrWeights[[i]]),nObs*Dim2Length*Dim3Length)))
            } 
          } else {
            ### check if 2D:
            if (length(dim(attrWeights[[i]]))==2) {
              ### 2D
              if(dim(attrWeights[[i]])[1]==nObs){
                attrWeightsM<-c(attrWeightsM,rep(t(attrWeights[[i]]),each=Dim3Length)) 
              } else {
                attrWeightsM<-c(attrWeightsM,rep(rep(attrWeights[[i]],nObs),each=Dim3Length))
              }
            } else {
              ### 3D
              if(dim(attrWeights[[i]])[1]==nObs){
                ### may or may not have 2nd dim
                if(dim(attrWeights[[i]])[2]==Dim2Length){
                  attrWeightsM<-c(attrWeightsM,c(aperm(attrWeights[[i]])))
                } else {
                  attrWeightsM<-c(attrWeightsM,c(t(matrix(rep(attrWeights[[i]],each=Dim2Length),ncol=Dim3Length))))
                }
              } else {
                if(dim(attrWeights[[i]])[2]==Dim2Length){
                  attrWeightsM<-c(attrWeightsM,rep(c(aperm(attrWeights[[i]])),nObs))
                } else{
                  attrWeightsM<-c(attrWeightsM,rep(attrWeights[[i]],nObs*Dim2Length))
                }
              }
            }
          }
        }
        attrWeightsM<-array(attrWeightsM,c(nObs*Dim2Length*Dim3Length,nAttrs))
      } else {
        attrWeightsM<-array(1/nAlts,c(nObs*Dim2Length*Dim3Length,nAttrs))
      }
      
      #### attrScalings
      ## could not be provided (=1)
      ## then for each attribute:
      ## could be attribute specific -> not mixed, some mixed, all mixed (2 or 3 dim)
      ## could be general -> not mixed, some mixed, all mixed (2 or 3 dim)
      ## 16 cases: attribute specific, choice specific, inter, intra...
      ## for each attribute:
      if(sum(lengths(attrScalings))!=1){
        attrScalingsM<-c()
        if(is.list(attrScalings)){
          ## attribute specific
          for(i in 1:nAttrs) {
            if(is.null(dim(attrScalings[[i]]))) {
              ### 1D, no mixing:
              if(length(attrScalings[[i]])==nObs){
                ### diff value for different obs
                attrScalingsM<-c(attrScalingsM,rep(c(attrScalings[[i]]),each=Dim2Length*Dim3Length))
              } else {
                attrScalingsM<-c(attrScalingsM,c(rep(c(attrScalings[[i]]),nObs*Dim2Length*Dim3Length)))
              } 
            } else {
              ### check if 2D:
              if (length(dim(attrScalings[[i]]))==2) {
                ### 2D
                if(dim(attrScalings[[i]])[1]==nObs){
                  attrScalingsM<-c(attrScalingsM,rep(t(attrScalings[[i]]),each=Dim3Length)) 
                } else {
                  attrScalingsM<-c(attrScalingsM,rep(rep(attrScalings[[i]],nObs),each=Dim3Length))
                }
              } else {
                ### 3D
                if(dim(attrScalings[[i]])[1]==nObs){
                  ### may or may not have 2nd dim
                  if(dim(attrScalings[[i]])[2]==Dim2Length){
                    attrScalingsM<-c(attrScalingsM,c(aperm(attrScalings[[i]])))
                  } else {
                    attrScalingsM<-c(attrScalingsM,c(t(matrix(rep(attrScalings[[i]],each=Dim2Length),ncol=Dim3Length))))
                  }
                } else {
                  if(dim(attrScalings[[i]])[2]==Dim2Length){
                    attrScalingsM<-c(attrScalingsM,rep(c(aperm(attrScalings[[i]])),nObs))
                  } else{
                    attrScalingsM<-c(attrScalingsM,rep(attrScalings[[i]],nObs*Dim2Length))
                  }
                }
              }
            }
          }
        } else {
          ## not attribute specific
          if(is.null(dim(attrScalings))) {
            ### 1D, no mixing:
            if(length(attrScalings)==nObs){
              ### diff value for different obs
              attrScalingsM<-c(attrScalingsM,rep(rep(c(attrScalings),each=Dim2Length*Dim3Length),nAlts))
            } else {
              attrScalingsM<-c(attrScalingsM,rep(c(rep(c(attrScalings),nObs*Dim2Length*Dim3Length)),nAlts))
            } 
          } else {
            ### check if 2D:
            if (length(dim(attrScalings))==2) {
              ### 2D
              if(dim(attrScalings)[1]==nObs){
                attrScalingsM<-c(attrScalingsM,rep(rep(t(attrScalings),each=Dim3Length),nAlts)) 
              } else {
                attrScalingsM<-c(attrScalingsM,rep(rep(rep(attrScalings,nObs),each=Dim3Length),nAlts))
              }
            } else {
              ### 3D
              if(dim(attrScalings)[1]==nObs){
                ### may or may not have 2nd dim
                if(dim(attrScalings)[2]==Dim2Length){
                  attrScalingsM<-c(attrScalingsM,rep(c(aperm(attrScalings)),nAlts))
                } else {
                  attrScalingsM<-c(attrScalingsM,c(rep(t(matrix(rep(attrScalings,each=Dim2Length),ncol=Dim3Length))),nAlts))
                }
              } else {
                if(dim(attrScalings)[2]==Dim2Length){
                  attrScalingsM<-c(attrScalingsM,rep(rep(c(aperm(attrScalings)),nObs),nAlts))
                } else{
                  attrScalingsM<-c(attrScalingsM,rep(rep(attrScalings,nObs*Dim2Length),nAlts))
                }
              }
            }
          }
        }
        ### build matrix
        attrScalingsM<-array(attrScalingsM,c(nObs*Dim2Length*Dim3Length,nAttrs*nAlts))
      } else {
        attrScalingsM<-array(1,c(nObs*Dim2Length,nAttrs))
        if (attrScalings!=1) warning("If you are not using attrScalings, please set it to 1")
      }
      
      ### process parameters:
      ### 8 cases:
      ### 3rd, 2nd, 1st dims different
      if(is.null(dim(erv))){
        if(length(erv)==nObs) {
          erv = rep(erv,each=Dim2Length*Dim3Length)
        } else {
          erv = rep(erv,Dim2Length*Dim3Length*nObs)
        }
      } else {
        if(length(dim(erv))==2){
          if(dim(erv)[1]==nObs){
            erv = rep(t(erv),each=Dim3Length)
          } else {
            erv = rep(rep(erv,each=Dim3Length),nObs) #doesn't matter here which rep first...
          }
        } else {
          ## dim ==3, 4 more cases
          if(dim(erv)[2]==Dim2Length){
            if(dim(erv)[1]==nObs){
              erv<-c(aperm(erv))
            } else {
              erv<-rep(c(aperm(erv)),nObs)
            }
          } else {
            if(dim(erv)[1]==nObs){
              erv<-c(t(matrix(rep(erv,each=Dim2Length),ncol=Dim3Length)))
            } else {
              erv<-rep(erv,nObs*Dim2Length)
            }
          }
        } 
      }
      
      if(is.null(dim(ts))){
        if(length(ts)==nObs) {
          ts = rep(ts,each=Dim2Length*Dim3Length)
        } else {
          ts = rep(ts,Dim2Length*Dim3Length*nObs)
        }
      } else {
        if(length(dim(ts))==2){
          if(dim(ts)[1]==nObs){
            ts = rep(t(ts),each=Dim3Length)
          } else {
            ts = rep(rep(ts,each=Dim3Length),nObs) #doesn't matter here which rep first...
          }
        } else {
          ## dim ==3, 4 more cases
          if(dim(ts)[2]==Dim2Length){
            if(dim(ts)[1]==nObs){
              ts<-c(aperm(ts))
            } else {
              ts<-rep(c(aperm(ts)),nObs)
            }
          } else {
            if(dim(ts)[1]==nObs){
              ts<-c(t(matrix(rep(ts,each=Dim2Length),ncol=Dim3Length)))
            } else {
              ts<-rep(ts,nObs*Dim2Length)
            }
          }
        } 
      }
      
      if(is.null(dim(phi1))){
        if(length(phi1)==nObs) {
          phi1 = rep(phi1,each=Dim2Length*Dim3Length)
        } else {
          phi1 = rep(phi1,Dim2Length*Dim3Length*nObs)
        }
      } else {
        if(length(dim(phi1))==2){
          if(dim(phi1)[1]==nObs){
            phi1 = rep(t(phi1),each=Dim3Length)
          } else {
            phi1 = rep(rep(phi1,each=Dim3Length),nObs) #doesn't matter here which rep first...
          }
        } else {
          ## dim ==3, 4 more cases
          if(dim(phi1)[2]==Dim2Length){
            if(dim(phi1)[1]==nObs){
              phi1<-c(aperm(phi1))
            } else {
              phi1<-rep(c(aperm(phi1)),nObs)
            }
          } else {
            if(dim(phi1)[1]==nObs){
              phi1<-c(t(matrix(rep(phi1,each=Dim2Length),ncol=Dim3Length)))
            } else {
              phi1<-rep(phi1,nObs*Dim2Length)
            }
          }
        } 
      }
      
      if(is.null(dim(phi2))){
        if(length(phi2)==nObs) {
          phi2 = rep(phi2,each=Dim2Length*Dim3Length)
        } else {
          phi2 = rep(phi2,Dim2Length*Dim3Length*nObs)
        }
      } else {
        if(length(dim(phi2))==2){
          if(dim(phi2)[1]==nObs){
            phi2 = rep(t(phi2),each=Dim3Length)
          } else {
            phi2 = rep(rep(phi2,each=Dim3Length),nObs) #doesn't matter here which rep first...
          }
        } else {
          ## dim ==3, 4 more cases
          if(dim(phi2)[2]==Dim2Length){
            if(dim(phi2)[1]==nObs){
              phi2<-c(aperm(phi2))
            } else {
              phi2<-rep(c(aperm(phi2)),nObs)
            }
          } else {
            if(dim(phi2)[1]==nObs){
              phi2<-c(t(matrix(rep(phi2,each=Dim2Length),ncol=Dim3Length)))
            } else {
              phi2<-rep(phi2,nObs*Dim2Length)
            }
          }
        } 
      }
      
      choiceVar = rep(choiceVar,each=Dim2Length*Dim3Length)
      
    }
    
    
    
    ##### Get probabilities under DFT model
    
    if((functionality=="prediction")|(functionality=="raw")){
      P=list()
      for (j in 1:nAlts){
        if(Dims==1){
          P[[altnames[j]]] = sapply(1:nObs,function(i) calculateDFTProbs(j, c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))
        }
        
        if(Dims==2){
          P[[altnames[j]]] = matrix(c(sapply(1:(nObs*Dim2Length),function(i) calculateDFTProbs(j, c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))),nObs,Dim2Length,byrow=TRUE)
        }
        
        if(Dims==3){
          P[[altnames[j]]] = aperm(array(c(sapply(1:(nObs*Dim2Length*Dim3Length),function(i) calculateDFTProbs(j, c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))),c(Dim3Length,Dim2Length,nObs)))
        }
      }
      
      P <- lapply(P, function(p) {
        # change likelihood of excluded columns
        if(is.vector(p)) p[!rows]  <- ifelse(functionality=="prediction",NA,1)
        if(is.matrix(p)) p[!rows,] <- ifelse(functionality=="prediction",NA,1)
        if(is.array(p) & length(dim(p))==3) p[!rows,,] <- ifelse(functionality=="prediction",NA,1)
        return(p)
      })
      
      if(!choiceNA) {
        P[["chosen"]]= P[[altnames[[1]]]]*(alternatives[[altnames[[1]]]]==choiceVar)
        for (i in 2:nAlts) P[["chosen"]]= P[["chosen"]] + P[[altnames[[i]]]]*(alternatives[[altnames[[i]]]]==choiceVar)
      }
      
    } else {
      
      if(Dims==1){
        P = sapply(1:nObs,function(i) calculateDFTProbs(c(choiceVar[i]), c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))
        P[!rows]  <- 1
      }
      
      if(Dims==2){
        P = matrix(c(sapply(1:(nObs*Dim2Length),function(i) calculateDFTProbs(c(choiceVar[i]), c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))),nObs,Dim2Length,byrow=TRUE)
        P[!rows,] <- 1
      }
      
      if(Dims==3){
        P= aperm(array(c(sapply(1:(nObs*Dim2Length*Dim3Length),function(i) calculateDFTProbs(c(choiceVar[i]), c(attrValuesM[i,]), c(availM[i,]), c(altStartM[i,]), c(attrWeightsM[i,]), c(attrScalingsM[i,]),erv[i],ts[i],phi1[i],phi2[i],nAlts,nAttrs))),c(Dim3Length,Dim2Length,nObs)))
        P[!rows,,] <- 1
      }
    }
    
    if(functionality=="output"){
      
      # turn scalar availabilities into vectors
      for(i in 1:length(avail)) if(length(avail[[i]])==1) avail[[i]] <- rep(avail[[i]], nObs)
      
      # Construct summary table of availabilities and market share
      availprint = colSums(rows*matrix(unlist(avail), ncol = length(avail))) # number of times each alt is available
      choicematrix = matrix(0,nrow=4,ncol=length(altnames))
      choicematrix[1,] = availprint
      j=1
      while(j<= length(altnames)){
        choicematrix[2,j]=sum(choiceVar==altcodes[j] & rows) # number of times each alt is chosen
        j=j+1
      }
      choicematrix[3,] = choicematrix[2,]/sum(rows)*100 # market share
      choicematrix[4,] = choicematrix[2,]/choicematrix[1,]*100 # market share controlled by availability
      rownames(choicematrix) = c("Times available","Times chosen","Percentage of choice overall","Percentage of choice when available")
      colnames(choicematrix) = altnames
      
      
      
      ## write diagnostics to a file named "modelName_tempOutput.txt" in a temporary directory.
      #apollo_control <- tryCatch( get("apollo_inputs", parent.frame(), inherits=FALSE )$apollo_control,
      #                            error=function(e){
      #                              cat("apollo_dft could not retrieve apollo_control. No diagnostics in output.\n")
      #                              return(NA)
      #                            } )
      #if(!(length(apollo_control)==1 && is.na(apollo_control))){
      #  fileName <- paste(apollo_control$modelName, "_tempOutput.txt", sep="")
      #  fileName <- file.path(tempdir(),fileName)
      #  fileConn <- tryCatch( file(fileName, open="at"),
      #                        error=function(e){
      #                          cat('apollo_dft could not write diagnostics to temporary file. No diagnostics in output.\n')
      #                          return(NA)
      #                        })
      #  if(!anyNA(fileConn)){
      #    sink(fileConn)
      #    on.exit({if(sink.number()>0) sink(); close(fileConn)})
      #    if(apollo_control$noDiagnostics==FALSE){
      #      cat('Overview of choices for DFT model component:\n')
      #      print(round(choicematrix,2))
      #      cat("\n")}
      #    if(sum(choicematrix[4,]==0)>0) cat("Warning: some alternatives are never chosen in your data!\n")
      #    if(any(choicematrix[4,]==1)) cat("Warning: some alternatives are always chosen when available!\n")
      #  }
      #}
      
      content <- list(round(choicematrix,2))
      if(any(choicematrix[4,]==0)) content[[length(content) + 1]] <- "Warning: some alternatives are never chosen in your data!"
      if(any(choicematrix[4,]==1)) content[[length(content) + 1]] <- "Warning: some alternatives are always chosen when available!"
      if(avail_set==TRUE) content[[length(content)+1]] <- paste0("Warning: Availability not provided (or some elements are NA).",
                                                                 "\n", paste0(rep(" ",9),collapse=""),"Full availability assumed.")
      if(warn1) content[[length(content)+1]] <- paste0("Notice: Not all of the attributes given in \"attrValues\" are \n",
                                                       paste0(rep(" ",8), collapse=""), "named in \"attrScalings\" or \"attrWeights\". These will\n", 
                                                       paste0(rep(" ",8), collapse=""), "consequently be ignored.")
      if(warn2) content[[length(content)+1]] <- paste0("Notice: Not all of the alternatives given in \"altStart\" are\n",
                                                       paste0(rep(" ",8), collapse=""),"named in \"alternatives\". These will consequently be ignored.")
      if(warn3) content[[length(content)+1]] <- paste0("Notice: A list was not supplied for \"altStart\".\n",
                                                       paste0(rep(" ",8), collapse=""),"Starting values for all alternatives will be set to zero.")
      apolloLog <- tryCatch(get("apollo_inputs", parent.frame(), inherits=TRUE )$apolloLog, error=function(e) return(NA))
      apollo_addLog("Overview of choices for DFT model component:", content, apolloLog)
    }
    
    return(P)
    
  }
  
  
}


calculateDFTProbs<-function(choiceVar,attribs, avail, altStart, attrWeights, attrScalings , erv, ts, phi1, phi2, nAlts, nAttrs){
  
  M<-matrix(c(attribs),nAlts,nAttrs,byrow=TRUE)*attrScalings
  
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
  kp=(avail)*c(1:nAlts)
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

