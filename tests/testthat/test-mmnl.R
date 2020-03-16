context("test-mmnl")

test_that("mmnl works", {
  rm(list = ls())
  apollo_initialise()

  apollo_control = list(
    modelName ="MMNL",# Make sure to use a new name for every model
    modelDescr ="Mixed logit model with inter and intra draws",
    indivID   ="ID",  # Name of column in the database with each individual's ID
    mixing    = TRUE, # TRUE for models that include random parameters
    nCores    = 2
  )
  data("apollo_modeChoiceData", package="apollo")
  database = apollo_modeChoiceData
  database = subset(database,database$SP==1 & database$ID<=100)

  apollo_beta=c(asc_car  = 0.0000,
                asc_bus  =-2.5111,
                asc_air  = 0.1133,
                asc_rail =-0.3724,
                b_tt_mu  =-0.0144,
                b_tt_sgBe= 0.0059,
                b_tt_sgWi=-0.0032,
                b_access =-0.0271,
                b_cost   =-0.0673)
  apollo_fixed = c("asc_car")

  apollo_draws = list(
    interDrawsType = "halton",
    interNDraws    = 10,
    interUnifDraws = c(),
    interNormDraws = c("draws_tt_between"),
    intraDrawsType = "mlhs",
    intraNDraws    = 10,
    intraUnifDraws = c(),
    intraNormDraws = c("draws_tt_within")
  )
  apollo_randCoeff = function(apollo_beta, apollo_inputs){
    randcoeff = list()
    randcoeff[["b_tt"]] = b_tt_mu + b_tt_sgBe*draws_tt_between + b_tt_sgWi*draws_tt_within
    return(randcoeff)
  }

  apollo_inputs <- apollo_validateInputs(apollo_beta, apollo_fixed, database,
                                 apollo_control, apollo_HB=NA, apollo_draws,
                                 apollo_randCoeff, apollo_lcPars=NA)

  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P = list()
    V = list()
    V[['car']]  = asc_car  + b_tt*time_car                         + b_cost*cost_car
    V[['bus']]  = asc_bus  + b_tt*time_bus  + b_access*access_bus  + b_cost*cost_bus
    V[['air']]  = asc_air  + b_tt*time_air  + b_access*access_air  + b_cost*cost_air
    V[['rail']] = asc_rail + b_tt*time_rail + b_access*access_rail + b_cost*cost_rail
    mnl_settings = list(
      alternatives  = c(car=1, bus=2, air=3, rail=4),
      avail         = list(car=av_car, bus=av_bus, air=av_air, rail=av_rail),
      choiceVar     = choice,
      V             = V
    )
    P[['model']] = apollo_mnl(mnl_settings, functionality)
    P = apollo_avgIntraDraws(P, apollo_inputs, functionality)
    P = apollo_panelProd(P, apollo_inputs, functionality)
    P = apollo_avgInterDraws(P, apollo_inputs, functionality)
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }


  model = apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)
  apollo_modelOutput(model)
  apollo_saveOutput(model)

  expect_equal(model$maximum, -1066.478941, tolerance=0.1, scale=1)
})
