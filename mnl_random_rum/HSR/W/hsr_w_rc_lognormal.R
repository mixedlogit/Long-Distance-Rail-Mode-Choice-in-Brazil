# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\6_Mixed model - SP Rail\\Model")

### Clear memory
rm(list = ls())

### Load Apollo library
install.packages('apollo')

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

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'rc_rum 2',
  modelDescr   = "MMNL with lognormal distribution",
  indivID      = "ID",
  mixing       = TRUE,
  nCores       = 6)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_hsr    = 0,
                asc_air    = 0,
                mu_btt     = -3,
                sigma_btt  = 1,
                b_fr     = 0,
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
  interDrawsType = "mlhs",
  interNDraws    = 500,
  interNormDraws = c("draws_tt",'draws_fr')
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_tt"]] = exp(mu_btt + sigma_btt*draws_tt) 

  return(randcoeff)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

set.seed(42)
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

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

#likelihood ratio test
apollo_lrTest('mnl_rum 4','mnl_rum 3')


#VTTS RUM
apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="lognormal",parName1='mu_btt',parName2="sigma_btt"))
apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="lognormal",parName1='mu_bfr',parName2="sigma_bfr"))

apollo_deltaMethod(model = model_rum,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_co"))[1]*60

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
                                 searchStart_settings=list(nCandidates=10))



### Print outputs of additional diagnostics to new output file (remember to close file writing when complete)
sink(paste(model$apollo_control$modelName,"_additional_output.txt",sep=""),split=TRUE)

# ----------------------------------------------------------------- #
#---- CONDITIONALS AND UNCONDITIONALS                            ----
# ----------------------------------------------------------------- #


plot(density(as.vector(unconditionals[["b_tt"]])))
plot(density(as.vector(unconditionals[["b_fr"]])))

conditionals <- apollo_conditionals(model,apollo_probabilities, apollo_inputs)

plot(density(as.vector(conditionals[["b_fr"]])))

mean(unconditionals[["b_tt"]])

sd(unconditionals[["b_tt"]])

summary(conditionals[["b_tt"]])

income_n = apollo_firstRow(database$hh_inc_abs, apollo_inputs)

summary(lm(conditionals[["v_tt"]][,2]~income_n))

write.csv(conditionals,paste(model$apollo_control$modelName,"_conditionals.csv",sep=""))

########## VTT ##################
unconditionals <- apollo_unconditionals(model,apollo_probabilities, apollo_inputs)
unconditionals$b_tt

#unconditionals
b_cost = model$estimate[['b_co']]
unconditionals$b_tt

#simulation
N = 1000
b_tt = -exp(rnorm(N, mean=model$estimate["mu_btt"], sd=model$estimate["sigma_btt"]))
b_cost = model$estimate[['b_co']]
wtp = b_tt/b_cost*60

mean(wtp)
median(wtp)
sd(wtp)

apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="ratio",parName1='b_fr',parName2="b_co"))

summary(wtp)

ggplot()+
  geom_density(wt)

plot(x = wtp)


a <- model$hessian
b <- model$varcovBGW
c <- model$varcov

a[c('mu_btt','sigma_btt','b_co'),c('mu_btt','sigma_btt','b_co')]
b[c('mu_btt','sigma_btt','b_co'),c('mu_btt','sigma_btt','b_co')]
c[c('mu_btt','sigma_btt','b_co'),c('mu_btt','sigma_btt','b_co')]

a
