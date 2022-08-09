#' Converts data from long to wide format.
#'
#' Converts choice data from long to wide format, with one row per observation as opposed to one row per alternative/observation.
#'
#' @param longData data.frame. Data in long format.
#' @param longToWide_settings List. Contains settings for this function. User input is required for all settings. 
#'                                \itemize{
#'                                  \item \strong{code{alternative_column}}: Character. Name of column in long data that contains the names of the alternatives (either numeric or character).
#'                                  \item \strong{code{alternative_specific_attributes}}: Character vector. Names of columns in long data with attributes that vary across alternatives within an observation.
#'                                  \item \strong{code{choice_column}}: Character. Name of column in long data that contains the choice.
#'                                  \item \strong{code{ID_column}}: Character. Name of column in long data that contains the ID of individuals.
#'                                  \item \strong{code{observation_column}}: Character. Name of column in long data that contains the observation index.
#'                                }
#' @return Silently returns a data.frame with the wide format version of the data.
#'         An overview of the data is printed to screen.
#' @export
apollo_longToWide=function(longData,longToWide_settings){
  ### CHECKS
  if(!is.data.frame(longData)) stop("The data object passed to \'apollo_longToWide\' needs to be a data.frame!")
  if(!is.list(longToWide_settings)) stop("The \'longToWide_settings\' argument passed to \'apollo_longToWide\' needs to be a list!")
  mandatory <- c("alternative_column", "alternative_specific_attributes", "choice_column", "ID_column", "observation_column")
  for(i in mandatory) if(!(i %in% names(longToWide_settings))) stop('The \'longToWide_settings\' list needs to include an object called "', i,'"!')
  for(i in names(longToWide_settings)[names(longToWide_settings)!="alternative_specific_attributes"]) if(!longToWide_settings[[i]]%in%colnames(longData)) stop('No column entitled "', longToWide_settings[[i]],'" found in data!')
  for(i in longToWide_settings$alternative_specific_attributes) if(!i%in%colnames(longData)) stop('No column entitled "', i,'" found in data!')
  if(length(longToWide_settings$alternative_specific_attributes)!=length(unique(longToWide_settings$alternative_specific_attributes))) stop("Some of the attribute names are duplicated in \'alternative_specific_attributes\'.")
  ### working copies
  alternative_column=longToWide_settings$alternative_column
  alternative_specific_attributes=longToWide_settings$alternative_specific_attributes
  choice_column=longToWide_settings$choice_column
  ID_column=longToWide_settings$ID_column
  observation_column=longToWide_settings$observation_column
  ### DATA SUMMARY
  alt=longData[,alternative_column]
  choice=longData[,choice_column]
  #overview=matrix(0,nrow=max(alt),ncol=2)
  overview=matrix(0,nrow=length(unique(alt)),ncol=2)
  rownames(overview)=paste0("Alternative ",unique(alt))
  colnames(overview)=c("Available","Chosen")
  for(j in 1:length(unique(longData$alt))){
    overview[j,1]=sum(longData$alt==unique(longData$alt)[j])
    #overview[j,2]=sum(alt*choice==j)
    overview[j,2]=sum(subset(longData,longData$alt==unique(longData$alt)[j])$choice)
  }
  ### GENERATED NEW COLUMNS
  #create ID-observation index
  longData$ID_obs=longData[,ID_column]+longData[,observation_column]/(1+max(longData[,observation_column]))
  ### further checks
  test1=length(unique(longData$ID_obs*choice))==(length(unique(longData$ID_obs))+1)
  test2=sum(overview[,2])==(length(unique(longData$ID_obs)))
  if(!(test1&test2)) stop("A single alternative needs to be chosen in each separate choice scenario!")
  ###sort by ID_obs
  longData=longData[order(longData$ID_obs),]
  #names of alts
  alternative_names=unique(longData[,alternative_column])
  #generic columns
  generic_columns=colnames(longData)[!colnames(longData)%in%c(alternative_specific_attributes,alternative_column,choice_column)]
  #take a subset of the data that uses only the generic part, and only once per observation
  generic_part=longData[,generic_columns]
  generic_part=subset(generic_part,sequence(rle(generic_part$ID_obs)$lengths)==1)
  ### alt specific
  alt_specific_part=longData[,c("ID_obs",alternative_specific_attributes,alternative_column,choice_column)]
  ### CREATE WIDE DATA
  wide=generic_part
  for(j in alternative_names){
    alt_part=subset(alt_specific_part,alt_specific_part$alt==j)
    colnames(alt_part)[2:ncol(alt_part)]=paste0(colnames(alt_part[2:ncol(alt_part)]),"_",j)
    alt_part[,paste0("avail_",j)]=1
    ###check which rows are missing, i.e. which ID_obs are not in the unique ones from the generic part
    ###append rows with 0s for those, for all non-ID_obs columns
    if(any(!(generic_part$ID_obs)%in%alt_part$ID_obs)){
      missing=generic_part$ID_obs[!(generic_part$ID_obs%in%alt_part$ID_obs)]
      tmp=matrix(0,nrow=length(missing),ncol=ncol(alt_part))
      colnames(tmp)=colnames(alt_part)
      tmp[,1]=missing
      alt_part=rbind(alt_part,tmp)
    }
    ### reorder
    alt_part=alt_part[order(alt_part$ID_obs),]
    ### add to main data
    wide=cbind(wide,alt_part[,2:ncol(alt_part)])
  }
  ###remove alt columns
  wide=wide[,!(names(wide) %in% paste0(alternative_column,"_",alternative_names))]
  ###create new choice column
  #wide$choice_new=as.matrix(wide[,(names(wide) %in% paste0(choice_column,"_",alternative_names))])%*%as.matrix(alternative_names)
  wide$choice_new=as.matrix(wide[,(names(wide) %in% paste0(choice_column,"_",alternative_names))])%*%as.matrix(1:length(alternative_names))
  ###remove individual choice columns
  wide=wide[,!(names(wide) %in% paste0(choice_column,"_",alternative_names))]
  ###remove ID_obs column
  wide=wide[,!(names(wide)=="ID_obs")]
  ### RETURN WIDE FORMAT DATA
  cat("\nData successfully turned from long to wide format.")
  cat("\n\nLong data had",length(alternative_names),"unique alternatives")
  cat("\nAlternative names were:",alternative_names)
  cat("\n\nThe following notation was used for the wide format:")
  #cat("\nFor alternative j, the following naming convention has been adopted:")
  #cat("\n    Attributes:   ",paste0(longToWide_settings$alternative_specific_attributes,"_j"))
  #cat("\n    Availability:  avail_j")
  overview_new=matrix(0,nrow=length(alternative_names),ncol=length(alternative_specific_attributes)+2)
  rownames(overview_new)=paste0("Alternative ",alternative_names)
  colnames(overview_new)=c(alternative_specific_attributes,"Availability","choice_new")
  for(j in 1:length(alternative_names)){
    overview_new[j,1:length(alternative_specific_attributes)]=paste0(alternative_specific_attributes,"_",alternative_names[j])
    overview_new[j,length(alternative_specific_attributes)+1]=paste0("avail_",alternative_names[j])
    overview_new[j,length(alternative_specific_attributes)+2]=j
  }
  overview_new=as.data.frame(overview_new)
  cat("\nThe choice variable is included in column choice_new.")
  cat("\n")
  print(overview_new)
  cat("\n\nOverview of choice data\n")
  print(overview)
  return(wide)
}

