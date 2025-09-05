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
# autoregressive  <- FALSE

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
## Different reallocat probababilities
## according to baseline partition at time = 1 
## but constant over time

rm(list = ls())
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
autoregressive  <- FALSE
source("0_priorSimulation.R")


########################################
## Prior simulations using 
## latent autoregressive model
## for reallocation probabilities
########################################
## Simulation 1 
## stationary process phi in (-1, 1)
## phi = 0.5
##################
rm(list = ls())

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition

name <- "differentEta_phi05_AR"
FirstPart= rep(c(1,2,3,4), each = 5)
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)
## small variance (sd = 0.5)
phi <- 0.5
kappa <- 0.5

nsim <- 5000
clusterSpecific <- FALSE
autoregressive <- TRUE
source("0_priorSimulation.R")

rm(list = ls())

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition

name <- "differentEta_phi08_AR"
FirstPart= rep(c(1,2,3,4), each = 5)
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)
## small variance (sd = 0.5)
phi <- 0.8
kappa <- 0.5

nsim <- 5000
clusterSpecific <- FALSE
autoregressive <- TRUE
source("0_priorSimulation.R")


########################################################################
## Simulation 2
## phi = 1.5
##################
rm(list = ls())
# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition
name <- "differentEta_phi15_AR"
FirstPart= rep(c(1,2,3,4), each = 5)
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)

## small variance (sd = 0.5)
phi <- 1.5
kappa <- 0.5

nsim <- 5000
clusterSpecific <- FALSE
autoregressive <- TRUE

source("0_priorSimulation.R")



########################################################################
## Simulation 3
## phi = -1.5 - negative
##################
rm(list = ls())

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition
# FirstPart= seq(1:20)
name <- "differentEta_phineg15_AR"
FirstPart= rep(c(1,2,3,4), each = 5)
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)
## stationary process phi in (-1, 1)
## small variance (sd = 0.5)
phi <- -1.5
kappa <- 0.5

nsim <- 5000
clusterSpecific <- FALSE
autoregressive <- TRUE

source("0_priorSimulation.R")

