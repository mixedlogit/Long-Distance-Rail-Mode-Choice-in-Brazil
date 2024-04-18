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

database <- database %>% 
  filter(PURP == 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'rc_rum 2',
  modelDescr   = "MNL RUM mode in WTP space",
  indivID      = "ID",
  mixing       = T,
  nCores       = 6
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_cr     = 0,
                asc_air    = 0,
                mu_btt        = 1.5,
                sigma_btt        = 1,
                mu_bttbus        = 1.5,
                sigma_bttbus        = 1,
                b_co          = 0,
                b_fr          = 0,
                b_age2_car     = 0,
                b_age2         = 0,
                b_age3         = 0,
                b_inc2         = 0,
                b_inc3_car     = 0,
                b_inc3         = 0
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_cr")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interNormDraws = c("draws_tt","draws_ttbus")
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_tt"]] = exp(mu_btt + sigma_btt*draws_tt) 
  randcoeff[["b_ttbus"]] = exp(mu_bttbus + sigma_bttbus*draws_ttbus) 
  
  return(randcoeff)
}

##
setwd(paste0(getwd(),"\\First Paper\\mnl_rum\\CR\\W"))
      
apollo_beta <- apollo_readBeta(apollo_beta,apollo_fixed,
                               inputModelName = "mnl_rum 15")

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
  
  V[['car']]  = asc_car + b_co*(b_tt*TT_CAR + CO_CAR) + b_age2_car*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3_car*INC_3
  
  V[['bus']]  = asc_bus +  b_co*(b_ttbus*TT_BUS + FA_BUS + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  V[['cr']] = asc_cr +  b_co*(b_tt*TT_CR +FA_CR + b_fr*FR_CR) 
  
  V[['air']] = asc_air + b_co*(b_tt*TT_AIR + FA_AIR + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
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

setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model/First Paper/mnl_rc_rum/CR/W/wtpsp log/")

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

