# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\Mixed model - SP Rail\\Model\\First Paper")

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database <-  read.csv('database_calib.csv',header = T)

database <- database %>% 
  filter(PURP == 1)

database <- database %>% 
  filter(PURP != 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'rc_rum 15',
  modelDescr   = "MNL RUM with 1 normal random parameter",
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
                asc_cr     = 0,
                asc_hsr    = 0,
                asc_air    = 0,
                mu_b_frbusair    = 0,
                sigma_b_frbusair   = 0,
                mu_b_frcrhsr    = 0,
                sigma_b_frcrhsr   = 0
                )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_hsr")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "mlhs",
  interNDraws    = 250,
  interNormDraws = c('draws_frbusair','draws_frcrhsr')
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_frbusair"]] = mu_b_frbusair + sigma_b_frbusair*draws_frbusair
  randcoeff[["b_frcrhsr"]] = mu_b_frcrhsr + sigma_b_frcrhsr*draws_frcrhsr

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
  
  V[['car']]  = asc_car  
  
  V[['bus']]  = asc_bus + (b_frbusair*FR_BUS) 
  
  V[['cr']] = asc_cr + (b_frcrhsr*FR_CR)
  
  V[['hsr']] = asc_hsr + (b_frcrhsr*FR_HSR)
  
  V[['air']] = asc_air + (b_frbusair*FR_AIR) 
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(car=1,bus=2,cr=3,hsr=4,air=5),
    avail        = list(car=DRIV_LIC,bus=BUS_AV,cr=CR_AV,hsr=HSR_AV,air=AIR_AV),
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
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\Mixed model - SP Rail\\Model\\First Paper\\rc_rum\\W\\normal")

model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #
apollo_modelOutput(model, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

apollo_saveOutput(model,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

apollo_combineResults()
# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

#likelihood ratio test
maxLik::
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

