## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# ####################################################### #
#### 1. Definition of core settings                        
# ####################################################### #

### Clear memory
rm(list = ls())

### Load libraries
library(apollo)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName  ="MNL",
  modelDescr ="Simple MNL model on mode choice SP data",
  indivID    ="ID"
)

# ####################################################### #
#### 2. Data loading                                   ####
# ####################################################### #

data("apollo_modeChoiceData", package="apollo")
database = apollo_modeChoiceData
rm(apollo_modeChoiceData)

### Use only SP data
database = subset(database,database$SP==1)

# ####################################################### #
#### 3. Parameter definition                           ####
# ####################################################### #

### Vector of parameters, including any that are kept fixed 
### during estimation
apollo_beta=c(asc_car  = 0,
              asc_bus  = 0,
              asc_air  = 0,
              asc_rail = 0,
              b_tt_car = 0,
              b_tt_bus = 0,
              b_tt_air = 0,
              b_tt_rail= 0,
              b_c      = 0)

### Vector with names (in quotes) of parameters to be
###  kept fixed at their starting value in apollo_beta.
### Use apollo_beta_fixed = c() for no fixed parameters.
apollo_fixed = c("asc_car")

# ####################################################### #
#### 4. Input validation                               ####
# ####################################################### #

apollo_inputs = apollo_validateInputs()

# ####################################################### #
#### 5. Likelihood definition                          ####
# ####################################################### #

apollo_probabilities=function(apollo_beta, apollo_inputs, 
                              functionality="estimate"){

  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as
  ### in mnl_settings, order is irrelevant.
  V = list()
  V[['car']] = asc_car + b_tt_car *time_car + b_c*cost_car
  V[['bus']] = asc_bus + b_tt_bus *time_bus + b_c*cost_bus
  V[['air']] = asc_air + b_tt_air *time_air + b_c*cost_air
  V[['rail']]= asc_rail+ b_tt_rail*time_rail+ b_c*cost_rail
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(car=1, bus=2, air=3, rail=4), 
    avail         = list(car=av_car, bus=av_bus, 
                         air=av_air, rail=av_rail), 
    choiceVar     = choice,
    V             = V
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)

  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)

  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ####################################################### #
#### 6. Model estimation and reporting                 ####
# ####################################################### #

model = apollo_estimate(apollo_beta, apollo_fixed, 
                        apollo_probabilities, 
                        apollo_inputs)

apollo_modelOutput(model)

#apollo_saveOutput(model)

# ####################################################### #
#### 7. Postprocessing of results                      ####
# ####################################################### #

### Use the estimated model to make predictions
predictions_base = apollo_prediction(model, 
                                     apollo_probabilities, 
                                     apollo_inputs)

### Now imagine the cost for rail increases by 10% 
### and predict again
database$cost_rail = 1.1*database$cost_rail
apollo_inputs   = apollo_validateInputs()
predictions_new = apollo_prediction(model, 
                                    apollo_probabilities, 
                                    apollo_inputs)

### Compare predictions
change=(predictions_new-predictions_base)/predictions_base
### Not interested in chosen alternative now, 
### so drop last column
change=change[,-ncol(change)]
### Summary of changes (possible presence of NAs due to
### unavailable alternatives)
summary(change)


## ---- eval=FALSE--------------------------------------------------------------
#  # ####################################################### #
#  #### 1. Definition of core settings
#  # ####################################################### #
#  
#  ### Clear memory
#  rm(list = ls())
#  
#  ### Load libraries
#  library(apollo)
#  
#  ### Initialise code
#  apollo_initialise()
#  
#  ### Set core controls
#  apollo_control = list(
#    modelName  ="MMNL",
#    modelDescr ="Simple MMNL model on mode choice SP data",
#    indivID    ="ID",
#    mixing     = TRUE,
#    nCores     = 2
#  )
#  
#  # ####################################################### #
#  #### 2. Data loading                                   ####
#  # ####################################################### #
#  
#  data("apollo_modeChoiceData", package="apollo")
#  database = apollo_modeChoiceData
#  rm(apollo_modeChoiceData)
#  
#  ### Use only SP data
#  database = subset(database,database$SP==1)
#  
#  ### Create new variable with average income
#  database$mean_income = mean(database$income)
#  
#  # ####################################################### #
#  #### 3. Parameter definition                           ####
#  # ####################################################### #
#  
#  ### Vector of parameters, including any that are kept fixed
#  ### during estimation
#  apollo_beta=c(asc_car  = 0,
#                asc_bus  =-2,
#                asc_air  =-1,
#                asc_rail =-1,
#                mu_tt    =-4,
#                sigma_tt = 0,
#                b_c      = 0)
#  
#  ### Vector with names (in quotes) of parameters to be
#  ###  kept fixed at their starting value in apollo_beta.
#  ### Use apollo_beta_fixed = c() for no fixed parameters.
#  apollo_fixed = c("asc_car")
#  
#  ### Set parameters for generating draws
#  apollo_draws = list(
#    interDrawsType = "halton",
#    interNDraws    = 500,
#    interUnifDraws = c(),
#    interNormDraws = c("draws_tt"),
#    intraDrawsType = "halton",
#    intraNDraws    = 0,
#    intraUnifDraws = c(),
#    intraNormDraws = c()
#  )
#  
#  ### Create random parameters
#  apollo_randCoeff = function(apollo_beta, apollo_inputs){
#    randcoeff = list()
#  
#    randcoeff[["b_tt"]] = -exp(mu_tt + sigma_tt*draws_tt)
#  
#    return(randcoeff)
#  }
#  
#  # ####################################################### #
#  #### 4. Input validation                               ####
#  # ####################################################### #
#  
#  apollo_inputs = apollo_validateInputs()
#  
#  # ####################################################### #
#  #### 5. Likelihood definition                          ####
#  # ####################################################### #
#  
#  apollo_probabilities=function(apollo_beta, apollo_inputs,
#                                functionality="estimate"){
#  
#    ### Attach inputs and detach after function exit
#    apollo_attach(apollo_beta, apollo_inputs)
#    on.exit(apollo_detach(apollo_beta, apollo_inputs))
#  
#    ### Create list of probabilities P
#    P = list()
#  
#    ### List of utilities: these must use the same names as
#    ### in mnl_settings, order is irrelevant.
#    V = list()
#    V[['car']]  = asc_car  + b_tt*time_car  + b_c*cost_car
#    V[['bus']]  = asc_bus  + b_tt*time_bus  + b_c*cost_bus
#    V[['air']]  = asc_air  + b_tt*time_air  + b_c*cost_air
#    V[['rail']] = asc_rail + b_tt*time_rail + b_c*cost_rail
#  
#    ### Define settings for MNL model component
#    mnl_settings = list(
#      alternatives  = c(car=1, bus=2, air=3, rail=4),
#      avail         = list(car=av_car, bus=av_bus,
#                           air=av_air, rail=av_rail),
#      choiceVar     = choice,
#      V             = V
#    )
#  
#    ### Compute probabilities using MNL model
#    P[['model']] = apollo_mnl(mnl_settings, functionality)
#  
#    ### Take product across observation for same individual
#    P = apollo_panelProd(P, apollo_inputs, functionality)
#  
#    ### Average draws
#    P = apollo_avgInterDraws(P, apollo_inputs, functionality)
#  
#    ### Prepare and return outputs of function
#    P = apollo_prepareProb(P, apollo_inputs, functionality)
#    return(P)
#  }
#  
#  # ####################################################### #
#  #### 6. Model estimation and reporting                 ####
#  # ####################################################### #
#  
#  model = apollo_estimate(apollo_beta, apollo_fixed,
#                          apollo_probabilities,
#                          apollo_inputs)
#  
#  apollo_modelOutput(model)
#  
#  #apollo_saveOutput(model)
#  
#  # ####################################################### #
#  #### 7. Postprocessing of results                      ####
#  # ####################################################### #
#  
#  ### Use the estimated model to make predictions
#  predictions_base = apollo_prediction(model,
#                                       apollo_probabilities,
#                                       apollo_inputs)
#  
#  ### Now imagine the cost for rail increases by 10%
#  ### and predict again
#  database$cost_rail = 1.1*database$cost_rail
#  apollo_inputs   = apollo_validateInputs()
#  predictions_new = apollo_prediction(model,
#                                      apollo_probabilities,
#                                      apollo_inputs)
#  
#  ### Compare predictions
#  change=(predictions_new-predictions_base)/predictions_base
#  ### Not interested in chosen alternative now,
#  ### so drop last column
#  change=change[,-ncol(change)]
#  ### Summary of changes (possible presence of NAs due to
#  ### unavailable alternatives)
#  summary(change)

