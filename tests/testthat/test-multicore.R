context("test-multicore")

test_that("multicore estimation works", {
  rm(list = ls())
  apollo_initialise()
  apollo_control = list(
    modelName  = "MNL2Cores",
    modelDescr = "Simple MNL model on mode choice RP data",
    indivID    = "ID",
    nCores     = 2
  )

  data("apollo_modeChoiceData", package="apollo")
  database = apollo_modeChoiceData
  database = subset(database,database$RP==1)
  database$mean_income = mean(database$income)

  apollo_beta=c(asc_car   = 0,
                asc_bus   = 0,
                asc_air   = 0,
                asc_rail  = 0,
                b_tt_car  = 0,
                b_tt_bus  = 0,
                b_tt_air  = 0,
                b_tt_rail = 0,
                b_access  = 0,
                b_cost    = 0)
  apollo_fixed = c("asc_car")

  apollo_inputs <- apollo_validateInputs(apollo_beta, apollo_fixed, database,
                                 apollo_control, apollo_HB=NA, apollo_draws=NA,
                                 apollo_randCoeff=NA, apollo_lcPars=NA)

  choiceAnalysis_settings <- list(
    alternatives = c(car=1, bus=2, air=3, rail=4),
    avail        = list(car=database$av_car, bus=database$av_bus, air=database$av_air, rail=database$av_rail),
    choiceVar    = database$choice,
    explanators  = database[,c("female","business","income")]
  )

  apollo_choiceAnalysis(choiceAnalysis_settings, apollo_control, database)

  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P = list()
    V = list()
    V[['car']]  = asc_car  + b_tt_car  * time_car                        + b_cost * cost_car
    V[['bus']]  = asc_bus  + b_tt_bus  * time_bus  + b_access * access_bus  + b_cost * cost_bus
    V[['air']]  = asc_air  + b_tt_air  * time_air  + b_access * access_air  + b_cost * cost_air
    V[['rail']] = asc_rail + b_tt_rail * time_rail + b_access * access_rail + b_cost * cost_rail
    mnl_settings = list(
      alternatives  = c(car=1, bus=2, air=3, rail=4),
      avail         = list(car=av_car, bus=av_bus, air=av_air, rail=av_rail),
      choiceVar     = choice,
      V             = V
    )
    P[['model']] = apollo_mnl(mnl_settings, functionality)
    P = apollo_panelProd(P, apollo_inputs, functionality)
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }

  model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)


  apollo_modelOutput(model)
  apollo_saveOutput(model)

  expect_equal(model$maximum, -1025.756379 , tolerance=0.1, scale=1)
})
