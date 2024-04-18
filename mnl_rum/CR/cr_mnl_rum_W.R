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

database <-  read.csv('database.cr_calib.csv',header = T)
database <-  read.csv('database.cr_valid.csv',header = T)

database <- database %>% 
  filter(PURP == 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'mnl_rum 17',
  modelDescr   = "MNL RUM model for choice mode data",
  indivID      = "ID")

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_cr     = 0,
                asc_air    = 0,
                b_tt          = 0,
                b_co          = 0,
                b_cocr        = 0,
                b_coair       = 0,
                b_fr          = 0,
                b_age2         = 0,
                b_age2_car     = 0,
                b_age3         = 0,
                b_inc2         = 0,
                b_inc3         = 0,
                b_inc3_car     = 0
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

  V[['car']]  = asc_car + (b_tt*TT_CAR + b_co*CO_CAR) + b_age2_car*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3_car*INC_3
  
  V[['bus']]  = asc_bus +  (b_tt*TT_BUS + b_co*FA_BUS + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  V[['cr']] = asc_cr + (b_tt*TT_CR + b_cocr*FA_CR + b_fr*FR_CR) 
 
  V[['air']] = asc_air +  (b_tt*TT_AIR + b_coair*FA_AIR + b_fr*FR_AIR) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
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
apollo_lrTest('mnl_rum 17','mnl_rum 15')

###combinar resultados
apollo_combineResults()

model_rum <- read_rds("mnl_rum 15_model.rds")

#VTTS RUM
apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_tt',parName2="b_co", multPar1=60))
apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_tt',parName2="b_cocr", multPar1=60))
apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_tt',parName2="b_coair", multPar1=60))

apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_co"))
apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_cocr"))
apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_coair"))

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
                  sharesTest_settings = list(alternatives = c(car=1,bus=2,cr=3,air=5), 
                                             choiceVar=database$CHOICE),
                                             subsamples=list(inc2=(database$INC_2==1),+inc0=(database$INC_2==0)))

citation("overdisp")
