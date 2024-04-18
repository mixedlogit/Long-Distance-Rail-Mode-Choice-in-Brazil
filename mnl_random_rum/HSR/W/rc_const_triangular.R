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

#trabalho
database <- database %>% 
  filter(PURP == 1)

#lazer
database <- database %>% 
  filter(PURP != 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'rc_rum 6',
  modelDescr   = "MMNL model with 1 random parameter",
  indivID      = "ID",
  mixing       = TRUE,
  nCores       = 6
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_hsr    = 0,
                asc_air    = 0,
                bound_btt      = -0.05,
                bound_bfr           = 0,
                b_co           = 0,
                b_age2         = 0,
                b_age2_car     = 0,
                b_age3_car     = 0,
                b_age3_bus     = 0,
                b_age3_air     = 0,
                b_inc2     = 0,
                b_inc2_air     = 0,
                b_inc3     = 0,
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
  
  V[['car']]  = asc_car + (b_tt*TT_CAR + b_co*CO_CAR) + b_age2_car*AGE_2 + b_age3_car*AGE_3 + b_inc2*INC_2 +
    b_inc3_car*INC_3
  
  V[['bus']]  = asc_bus +  (b_tt*TT_BUS + b_co*FA_BUS + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3_bus*AGE_3 + 
    b_inc2*INC_2 + b_inc3_car*INC_3
  
  V[['hsr']] = asc_hsr + (b_tt*TT_HSR + b_co*FA_HSR + b_fr*FR_HSR)
  
  V[['air']] = asc_air + (b_tt*TT_AIR +  b_co*FA_AIR + b_fr*FR_AIR) + b_age2*AGE_2 + b_age3_air*AGE_3 + 
    b_inc2_air*INC_2 + b_inc3*INC_3
  
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

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #
apollo_modelOutput(model, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

apollo_saveOutput(model,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model\\First Paper\\mnl_rc_rum/HSR/W/const_triang/")

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

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

#unconditionals
unconditionals <- apollo_unconditionals(model,apollo_probabilities, apollo_inputs)
conditionals <- apollo_conditionals(model,apollo_probabilities, apollo_inputs)

plot(density(as.vector(unconditionals[["b_tt"]])))
plot(density(as.vector(unconditionals[["b_fr"]])))
plot(density(as.vector(conditionals$b_tt$post.mean)))
plot(density(as.vector(conditionals$b_fr$post.mean)))


mean(unconditionals[["b_tt"]])

max(unconditionals[["b_tt"]])

sd(unconditionals[["b_tt"]])

summary(unconditionals[['b_tt']])
summary(conditionals[['b_fr']])
summary(conditionals)
