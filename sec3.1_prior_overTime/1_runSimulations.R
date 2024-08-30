########################################
## Some scenario for prior simulations
########################################
##################
## Simulation 1 
##################
# ## replicating one example frome PQD paper
# N <- 20;
# ## DP concentration parameter
# M <- 1;
# ## time
# ntime <- 10;
# ## Centering partition
# # FirstPart= seq(1:20)
# FirstPart <- NULL
# clusterSpecific <- FALSE

# nsim <- 5000

# alpha <- 0.75
# name <- "sameAlpha_075"
# source("0_priorSimulationR.R")

# alpha <- 0.5
# name <- "sameAlpha_05"
# source("0_priorSimulationR.R")

# alpha <- 0.9
# name <- "sameAlpha_09"
# source("0_priorSimulationR.R")
##################
## Simulation 2  
##################
## Different pegging probababilities
## according to baseline partition at time = 1 
## but constant over time

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition
# FirstPart= seq(1:20)
name <- "diffAlpha_constantOverTime"
FirstPart= rep(c(1,2,3,4), each = 5)
alpha <- rep(c(0.25, 0.5, 0.75 ,0.95), each = 5)
nsim <- 5000
clusterSpecific <- FALSE
source("0_priorSimulation.R")

##################
## Simulation 3 
##################
## Different pegging probababilities
## based on the prior partition clustering
## still constant over time

# N <- 20;
# ## DP concentration parameter
# M <- 1;
# ## time
# ntime <- 10;
# ## Centering partition
# # FirstPart= seq(1:20)
# name <- "clusterSpecificAlpha_constantOverTime"
# FirstPart= rep(c(1,2,3,4), each = 5)
# alpha <- rep(c(0.25, 0.5, 0.75 ,0.95), each = 5)
# nsim <- 5000
# clusterSpecific <- TRUE
# source("0_priorSimulationR.R")

##################
## Simulation 4 - high alpha
##################
## Different pegging probababilities
## according to baseline partition at time = 1 
## but constant over time

# # n of units to partition
# N <- 20;
# ## DP concentration parameter
# M <- 1;
# ## time
# ntime <- 10;
# ## Centering partition
# # FirstPart= seq(1:20)
# name <- "diffAlpha_constantOverTime_high"
# FirstPart= rep(c(1,2,3,4), each = 5)
# alpha <- rep(c(0.6, 0.7, 0.8 ,0.9), each = 5)
# nsim <- 5000
# clusterSpecific <- FALSE
# source("0_priorSimulationR.R")



