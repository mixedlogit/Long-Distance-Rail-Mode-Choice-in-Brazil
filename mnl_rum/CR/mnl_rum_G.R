# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model\\First Paper")

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database <-  read.csv('database_calib.csv',header = T)
database <-  read.csv('database_calib.csv',header = T)

database <-  read.csv('database.cr_calib.csv',header = T)
database <-  read.csv('database.cr_valid.csv',header = T)

database <-  read.csv('database.hsr_calib.csv',header = T)
database <-  read.csv('database.hsr_valid.csv',header = T)

database <- database %>% 
  filter(PURP == 1)

database <- database %>% 
  filter(PURP != 1)

database <- database %>% 
  mutate(CO_CAR = TO_CAR + PE_CAR)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'mnl_rum 1',
  modelDescr   = "MNL RUM model for choice mode data",
  indivID      = "ID"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_cr     = 0,
                asc_air    = 0,
                b_ttcar        = 0,
                b_ttbus        = 0,
                b_ttcr         = 0,
                b_ttair        = 0,
                b_cocar        = 0,
                b_cobus        = 0,
                b_cocr         = 0,
                b_coair        = 0,
                b_fr           = 0
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_cr")

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

  V[['car']]  = asc_car + (b_ttcar*TT_CAR + b_cocar*CO_CAR)
  
  V[['bus']]  = asc_bus +  (b_ttbus*TT_BUS + b_cobus*FA_BUS + b_fr*FR_BUS) 
  
  V[['cr']] = asc_cr + (b_ttcr*TT_CR + b_cocr*FA_CR + b_fr*FR_CR) 
 
  V[['air']] = asc_air +  (b_ttair*TT_AIR + b_coair*FA_AIR + b_fr*FR_BUS)
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(car=1,bus=2,cr=3,air=5),
    avail        = list(car=DRIV_LIC,bus=BUS_AV,cr=CR_AV,air=AIR_AV),
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
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\SP Rail\\Rail Choice Data_Doutorado Cassiano\\First Paper\\mnl_rum\\W\\PFSP")

model_rum <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

apollo_modelOutput(model_rum, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

apollo_saveOutput(model_rum,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

#likelihood ratio test

apollo_lrTest('mnl_rum 24','mnl_rum 49')

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

