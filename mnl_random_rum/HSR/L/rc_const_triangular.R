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

database <- database %>% 
  filter(PURP != 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'rc_rum 2',
  modelDescr   = "MMNL model with random parameters",
  indivID      = "ID",
  mixing       = TRUE,
  nCores       = 5
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car     = 0,
                asc_bus     = 0,
                asc_hsr     = 0,
                asc_air     = 0,
                bound_btt      = 0,
                bound_bfr      = 0,
                b_cocar     = 0,
                b_co        = 0,
                b_coair     = 0,
                b_age2     = 0,
                b_age2_air     = 0,
                b_age3     = 0,
                b_inc2     = 0,
                b_inc3         = 0,
                b_inc3_car     = 0
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_hsr")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c("draws_tta",'draws_ttb',"draws_fra",'draws_frb')
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_tt"]] = bound_btt + bound_btt*(draws_tta + draws_ttb)
  randcoeff[["b_fr"]] = bound_bfr + bound_bfr*(draws_fra + draws_frb)
  
  return(randcoeff)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  
  V[['car']]  = asc_car + (b_tt*TT_CAR + b_cocar*CO_CAR) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3_car*INC_3
  
  V[['bus']]  = asc_bus +  (b_tt*TT_BUS + b_co*FA_BUS + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  V[['hsr']] = asc_hsr + (b_tt*TT_HSR + b_co*FA_HSR + b_fr*FR_HSR) 
  
  V[['air']] = asc_air +  (b_tt*TT_AIR + b_coair*FA_AIR + b_fr*FR_AIR) + b_age2_air*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
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

  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)

  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs)


model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs,
                        estimate_settings = list(constraints=("mu_btt-sigma_btt==0")))

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #
apollo_modelOutput(model, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

apollo_saveOutput(model,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model\\First Paper\\mnl_rc_rum/HSR/L/const_triang/")

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #
#unconditionals
unconditionals <- apollo_unconditionals(model,apollo_probabilities, apollo_inputs)
conditionals <- apollo_conditionals(model,apollo_probabilities, apollo_inputs)

plot(density(as.vector(unconditionals[["b_tt"]])))
plot(density(as.vector(unconditionals[["b_fr"]])))

plot(density(as.vector(conditionals$post.mean)))
plot(density(as.vector(conditionals$b_tt$post.mean)))
plot(density(as.vector(conditionals$b_fr$post.mean)))


mean(unconditionals[["b_tt"]])/model$estimate[['b_co']]*60

max(unconditionals[["b_tt"]])

sd(unconditionals[["b_tt"]])

summary(conditionals[['b_tt']])
summary(conditionals[['b_fr']])
summary(conditionals)

#likelihood ratio test
apollo_lrTest('mnl_rum 4','mnl_rum 3')

### Optional speedTest
#speedTest_settings=list(
#   nDrawsTry = c(50, 75, 100),
#   nCoresTry = 1:3,
#   nRep      = 10
#)

#apollo_speedTest(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, speedTest_settings)

## Optional: searching for starting value
apollo_beta = apollo_searchStart(apollo_beta,
                                 apollo_fixed,
                                 apollo_probabilities,
                                 apollo_inputs,
                                 searchStart_settings=list(nCandidates=5))


## VTT
vtt <- unconditionals$b_tt/model$estimate[['b_co']]
mean(vtt)*60
sd(vtt)

vtt <- unconditionals$b_tt/model$estimate[['b_cocar']]
mean(vtt)*60
sd(vtt)

vtt <- unconditionals$b_tt/model$estimate[['b_coair']]
mean(vtt)*60
sd(vtt)

vtt_fr <- unconditionals$b_fr/model$estimate[['b_co']]
mean(vtt_fr)*60
sd(vtt)

plot(density(vtt*60))

apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_co"))[1]*60
apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_coair"))[1]*60

