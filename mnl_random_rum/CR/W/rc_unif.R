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
  modelDescr   = "MMNL with all random parameters",
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
                mu_b_ttcb     = 0,
                sigma_b_ttcb  = 0,
                mu_b_ttcrhsr     = 0,
                sigma_b_ttcrhsr  = 0,
                mu_b_ttair     = 0,
                sigma_b_ttair  = 0,
                mu_b_fr     = 0,
                sigma_b_fr  = 0,
                b_cocb     = 0,
                b_cocrhsr  = 0,
                b_coair    = 0
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_hsr")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "mlhs",
  interNDraws    = 1000,
  interUnifDraws = c("draws_ttcb","draws_ttcrhsr","draws_ttair","draws_fr")
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_ttcb"]] = mu_b_ttcb + sigma_b_ttcb*(2*draws_ttcb-1) 
  randcoeff[["b_ttcrhsr"]] = mu_b_ttcrhsr + sigma_b_ttcrhsr*(2*draws_ttcrhsr-1) 
  randcoeff[["b_ttair"]] = mu_b_ttair + sigma_b_ttair*(2*draws_ttair-1) 
  randcoeff[["b_fr"]] = mu_b_fr + sigma_b_fr*(2*draws_fr-1) 
  
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
  
  V[['car']]  = asc_car + (b_ttcb*TT_CAR + b_cocb*CO_CAR)
  
  V[['bus']]  = asc_bus +  (b_ttcb*TT_BUS + b_cocb*FA_BUS + b_fr*FR_BUS) 
  
  V[['cr']] = asc_cr + (b_ttcrhsr*TT_CR + b_cocrhsr*FA_CR + b_fr*FR_CR)
  
  V[['hsr']] = asc_hsr + (b_ttcrhsr*TT_HSR + b_cocrhsr*FA_HSR + b_fr*FR_HSR)
  
  V[['air']] = asc_air + (b_ttair*TT_AIR + b_coair*FA_AIR + b_fr*FR_AIR) 
  
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
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\Mixed model - SP Rail\\Model\\First Paper\\rc_rum\\W\\uniforme")

model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #
apollo_modelOutput(model, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

apollo_saveOutput(model,saveOutput_settings = list(printPVal=T,printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

conditionals <- apollo_conditionals(model,apollo_probabilities,apollo_inputs)

ggplot()+
  geom_density(mapping = aes(post.mean),data = conditionals$b_fr)

unconditionals <- apollo_unconditionals(model,apollo_probabilities,apollo_inputs)

ggplot()+
  geom_density(mapping = aes(),data = asunconditionals$b_ttcb)

class(unconditionals$b_ttcb)
as.data.frame(unconditionals$b_ttcb)
a <- as.data.frame(unconditionals$b_ttcb)
mean(unconditionals$b_ttcb[2,])
sd(unconditionals$b_ttcb)
conditionals$b_ttcb$post.mean[2]

mean(conditionals$b_ttcb$post.mean)

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

