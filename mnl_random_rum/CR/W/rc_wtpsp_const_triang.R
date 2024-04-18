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
  modelName    = 'rc_rum 3',
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
                asc_cr     = 0,
                asc_air    = 0,
                bound_btt     = 1.5,
                b_co          = 0,
                b_cocr          = 0,
                b_coair          = 0,
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
  interUnifDraws = c("draws_tta",'draws_ttb')
)

### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["b_tt"]] = bound_btt + bound_btt*(draws_tta + draws_ttb)/2

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
  
  V[['car']]  = asc_car + b_co*(b_tt*TT_CAR + CO_CAR) + b_age2_car*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3_car*INC_3
  
  V[['bus']]  = asc_bus +  b_co*(b_tt*TT_BUS + FA_BUS + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  V[['cr']] = asc_cr +  b_cocr*(b_tt*TT_CR +FA_CR + b_fr*FR_CR) 
  
  V[['air']] = asc_air + b_coair*(b_tt*TT_AIR + FA_AIR + b_fr*FR_BUS) + b_age2*AGE_2 + b_age3*AGE_3 +
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

setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\4_Mixed model - SP Rail\\Model/First Paper/mnl_rc_rum/CR/W/wtpsp const triang/")

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

#likelihood ratio test
apollo_lrTest('rc_rum 2','rc_rum 3')

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

unconditionals <- apollo_unconditionals(model,apollo_probabilities, apollo_inputs)
conditionals <- apollo_conditionals(model,apollo_probabilities, apollo_inputs)

plot(density(as.vector(unconditionals[["b_tt"]])))
plot(density(as.vector(unconditionals[["b_fr"]])))

plot(density(as.vector(conditionals$post.mean)))

median(unconditionals[["b_tt"]])*60

sd(unconditionals[["b_tt"]])

summary(conditionals[["b_tt"]])

class(model$hessian[c("mu_btt","sigma_btt"),c("mu_btt","sigma_btt")])

m <- matrix(model$hessian[c("mu_btt","sigma_btt"),c("mu_btt","sigma_btt")],nrow=2)

m <- -1/m #negative inverse of the Fisher matrix
id <- diag(3)
id[1:2,1:2] <- m #matrix as in Bliemer and Rose (2013)
id

se <- vector(length = 10000) #number of draws
mean.b <- model$estimate[["mu_btt"]] #estimated mean of the WTP distribution (lognormal)
sg.b <- model$estimate[["sigma_btt"]] #estimated sd of the WTP distribution (lognormal)
z <- vector(length = length(se))
beta <- vector(length = length(se))

for (i in seq(se)) {
  z[i] <- rnorm(1,mean = 0,sd = 1)
  beta[i] <- exp(mean.b + sg.b*z[i])
  
  se[[i]] <- sqrt(c(beta[i],beta[i]*z[i],beta[i]*sg.b) %*% id %*% c(beta[i],beta[i]*z[i],beta[i]*sg.b))
}

mean(beta)
median(beta)
mean(beta)*60
median(beta)*60

mean.se <- mean(se)
ci <- c(mean(beta)-1.96*mean.se,mean(beta)+1.96*mean.se)
ci
plot(density(beta))

apollo_deltaMethod(model = model,deltaMethod_settings = list(operation="lognormal",parName1='mu_btt',parName2="sigma_btt"))


N = 1000000
#b_attr = rnorm(N, mean=model$estimate["mu_attr"], sd=model$estimate["sg_attr"])
#b_wtp = exp(rnorm(N, mean=model$estimate["mu_btt"], sd=model$estimate["sigma_btt"]))
#wtp_attr = b_attr/b_cost
mean(b_wtp)*60; sd(b_wtp)
median(b_wtp)*60

plot(density(b_wtp))

wtp_tt <- as.data.frame(model$estimate[["bound_btt"]] + model$estimate[["bound_btt"]]*(runif(N) + runif(N))/2)
wtp_ttbus <- as.data.frame(model$estimate[["bound_ttbus"]] + model$estimate[["bound_ttbus"]]*(runif(N) + runif(N))/2)
names(wtp_tt) <- "wtp_tt"
names(wtp_ttbus) <- "wtp_ttbus"


ggplot()+
  geom_density(mapping = aes(wtp_tt*60,col="TT"),data = wtp_tt)+
  geom_density(mapping = aes(wtp_ttbus*60,col="TT_BUS"),data = wtp_ttbus)+
  theme_minimal()

  plot(density(wtp_tt))
plot(density(wtp_ttbus))

mean(wtp_tt)
sd(wtp_tt)

write.csv(conditionals,paste(model$apollo_control$modelName,"_conditionals.csv",sep=""))

# ----------------------------------------------------------------- #
#---- switch off writing to file                                 ----
# ----------------------------------------------------------------- #

if(sink.number()>0) sink()