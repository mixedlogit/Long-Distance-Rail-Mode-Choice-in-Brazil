# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model")

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database <-  read.csv('database.hsr_calib.csv',header = T)
database <-  read.csv('database.hsr_valid.csv',header = T)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'mnl_rum 30',
  modelDescr   = "MNL RUM model for choice mode data",
  indivID      = "ID"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_hsr     = 0,
                asc_air     = 0,
                b_ttcb     = 0,
                b_tthsr     = 0,
                b_ttair     = 0,
                b_fr        = 0,
                b_cocb        = 0,
                b_cohsr       = 0,
                b_coair       = 0,
                b_age1car        = 0,
                b_age1        = 0,
                b_age2        = 0,
                b_age2bus        = 0,
                b_inc1car           = 0,
                b_inc1bus           = 0,
                b_inc1air           = 0,
                b_inc2car              = 0,
                b_inc2              = 0
                )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_hsr")

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()

  V[['car']]  = asc_car + (b_ttcb*TT_CAR + b_cocb*CO_CAR) + b_age1car*AGE_1 + b_age2*AGE_2 + b_inc1car*INC_1 + b_inc2car*INC_2
  
  V[['bus']]  = asc_bus +  (b_ttcb*TT_BUS + b_cocb*FA_BUS + b_fr*FR_BUS) + b_age1*AGE_1 + b_age2bus*AGE_2 + b_inc1bus*INC_1 + b_inc2*INC_2 
  
  V[['hsr']] = asc_hsr + (b_tthsr*TT_HSR + b_cohsr*FA_HSR + b_fr*FR_HSR) 
 
  V[['air']] = asc_air +  (b_ttair*TT_AIR + b_coair*FA_AIR + b_fr*FR_AIR) + b_age1*AGE_1 + b_age2*AGE_2 + b_inc1air*INC_1 + b_inc2*INC_2
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(car=1,bus=2,hsr=4,air=5),
    avail        = list(car=DRIV_LIC,bus=BUS_AV,hsr=HSR_AV,air=AIR_AV),
    choiceVar    = CHOICE,
    V            = V
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
  #defining panel structure
  P = apollo_panelProd(P,apollo_inputs,functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
 
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model_rum <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

apollo_modelOutput(model_rum, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))


apollo_saveOutput(model_rum,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model\\First Paper\\mnl_rum\\HSR\\G")


# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

#likelihood ratio test

apollo_lrTest('mnl_rum 30','mnl_rum 29')

##COMBINE RESULTS ####

apollo_combineResults(combineResults_settings = list(printClassical =TRUE, printPVal=TRUE))


#VTTS RUM
vtts_car <- model_rum$estimate[['b_tt']]/model_rum$estimate[['b_co']]*60
vtts_bus <- model_rum$estimate[['b_tt']]/model_rum$estimate[['b_co']]*60
vtts_air <- model_rum$estimate[['b_tt']]/model_rum$estimate[['b_co']]*60
vtts_rail <- model_rum$estimate[['b_tt']]/model_rum$estimate[['b_co']]*60

vtts_car
vtts_bus
vtts_air
vtts_rail

### Use the estimated model to make predictions
predictions_base = apollo_prediction(model_rum, apollo_probabilities, apollo_inputs)

#time
elast_ttcar <- (1-predictions_base[,'car'])*model_rum$estimate[['b_tt']]*database$time_car
agg_ttcar <- sum(elast_ttcar*(predictions_base[,'car']/sum(predictions_base[,'car'])))

elast_ttbus <- (1-predictions_base[,'bus'])*model_rum$estimate[['b_tt']]*database$time_bus
agg_ttbus <- sum(elast_ttbus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

elast_ttair <- (1-predictions_base[,'air'])*model_rum$estimate[['b_tt']]*database$time_air
agg_ttair <- sum(elast_ttair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

elast_ttrail <- (1-predictions_base[,'rail'])*model_rum$estimate[['b_tt']]*database$time_rail
agg_ttrail <- sum(elast_ttrail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))

#cost
elast_cocar <- (1-predictions_base[,'car'])*model_rum$estimate[['b_co']]*database$cost_car
agg_cocar <- sum(elast_cocar*(predictions_base[,'car']/sum(predictions_base[,'car'])))

elast_cobus <- (1-predictions_base[,'bus'])*model_rum$estimate[['b_co']]*database$cost_bus
agg_cobus <- sum(elast_cobus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

elast_coair <- (1-predictions_base[,'air'])*model_rum$estimate[['b_co']]*database$cost_air
agg_coair <- sum(elast_coair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

elast_corail <- (1-predictions_base[,'rail'])*model_rum$estimate[['b_co']]*database$cost_rail
agg_corail <- sum(elast_corail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))

#access
elast_accbus <- (1-predictions_base[,'bus'])*model_rum$estimate[['b_acc']]*database$access_bus
agg_accbus <- sum(elast_accbus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

elast_accair <- (1-predictions_base[,'air'])*model_rum$estimate[['b_acc']]*database$access_air
agg_accair <- sum(elast_accair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

elast_accrail <- (1-predictions_base[,'rail'])*model_rum$estimate[['b_acc']]*database$cost_rail
agg_accrail <- sum(elast_accrail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))


#enumerado
agg_ttcar
agg_ttbus
agg_ttair
agg_ttrail
agg_cocar
agg_cobus
agg_coair
agg_corail
agg_accbus
agg_accair
agg_accrail

### Use the estimated model to make predictions
predictions_base = apollo_prediction(model_rum, apollo_probabilities, apollo_inputs)

predictions_base <- as.data.frame(predictions_base)
predictions_base$log <- log(predictions_base$chosen)
sum(predictions_base$log)
write.csv(predictions_base,file = 'resultados_rum_alt.csv')

### Look at summary of the predicted choice probabilities
summary(predictions_base)

apollo_sharesTest(model_rum,apollo_probabilities,apollo_inputs,
                  sharesTest_settings = list(alternatives = c(auto=3,metro=1,bus=2), choiceVar=database$choice))

