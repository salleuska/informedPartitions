## This script sets up a simple simulation scenario to compare the 
## posterior results using different informative priors on the partition space 
## The model used is a mixture of normal distributions with different means and equal variance
## Here we use the centered Partition process in Paganin et al. and the informed partition
##############################
# see https://github.com/salleuska/CPLogit
library(CPLogit)  ## CP process model
library(salso)    ## to estimate the posterior partition
library(drpm)     ## informed partition - devtools::install_github("gpage2990/drpm")
library(mclust)   ## adjustedRandIndex
##############################
## script with functions that compute waic and lpml
source("utils/LSP_MH.R")
##############################
## script with functions that compute waic and lpml
source("utils/lpml.R")
##############################
## data simulation - replicating simulation study 2 in sec 4.1
## with the difference that the standard deviation is 1
## 4 clusters with 25 observation each
## different degrees of separation 
##############################
## Using slurm to run each scenario separately 
listScenari <- read.table(file = "simScenario.csv", header = T, sep = ",")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
###################################
## TMP - parameters to test script ##
# task_id <- 5
# niter <- 100
# nburn <- 50
# nthin <- 1
# nout <- (niter - nburn)/nthin
# index <- seq(nburn +1, niter, by = nthin)

# niterLSP <- 100
# burnLSP <- niterLSP/2
# keepLSP <- seq((burnLSP+1), niterLSP, by = 1)

# ndata <- 3
###################################
## MCMC SETTINGS - CP and IP model
niter <- 10000
nburn <- 5000
nthin <- 1
nout <- (niter - nburn)/nthin
index <- seq(nburn +1, niter, by = nthin)

## MCMC SETTINGS - LSP 
## since it is MH we use more iterations
niterLSP <- 50000
burnLSP <- niterLSP/2
keepLSP <- seq((burnLSP+1), niterLSP, by = 5)

#####################
scenario <- listScenari[task_id, ]
#####################
## Functions
###################
# Function to generate data 
time.dat.gen <- function(N, Tm, auto_cor=0, sd=1){
  ## @N  = number of observations
  ## @Tm = number of time points
  ## @auto_corr = correlation between time points
  ## @sd = standard deviation
  if(auto_cor == 0){
        YMat <- t(replicate(N, arima.sim(model=list(), n=Tm, sd=sd)))
  } else {
    YMat <- t(replicate(N, arima.sim(model=list(ar=c(auto_cor)), n=Tm, sd=sd)))
  }
  YMat
}
###################
## simulation setting
args <- c(100, NA, 100, 1, 0.0, 1)
ndata       <- as.numeric(args[1]) # number of replications (i.e. datasets to generate)
datatype    <- as.numeric(args[2]) # outdated
N           <- as.numeric(args[3]) # total number of observations
Tm          <- as.numeric(args[4]) # number of time points
correlation <- as.numeric(args[5]) # to be used with multiple time points
stdev       <- as.numeric(args[6])

cat("ndata = ", ndata, "\n")
cat("N = ", N, "\n")
cat("datatype = ", datatype, "\n")
cat("Tm = ", Tm, "\n")
cat("correlation = ", correlation, "\n")
cat("stdev = ", stdev, "\n")

##############################
set.seed(134);

## Set parameters for data simulation
## mean for each clusters
mn_sep <- scenario$meanSep
mn  <- seq(-mn_sep, 2*mn_sep, by=mn_sep)
## correlation between multiple time points for each cluster
cor <- rep(correlation, 4)
  
## true partition
trueP <- cbind(rep(1:4, each=N/4))

## set initial partition
if(scenario$initialPartition=="true") {
  initialP <- trueP # 
} else if(scenario$initialPartition=="merge") { 
  initialP <- rep(1:2, each=N/2) # 
} else if(scenario$initialPartition=="split") { 
  initialP <- c(rep(1, 12), rep(2, 13), rep(3, 12), rep(4, 13),
          rep(5, 12), rep(6, 13), rep(7, 12), rep(8, 13))
} else {stop("option for initial partition not supported")}

## data structures for results
## adjusted rand index between estimated partition and initial one
ariInitialPartCP <- matrix(NA, nrow=ndata, ncol=6)
ariInitialPartIP <- matrix(NA, nrow=ndata, ncol=6)
ariInitialPartLSP <- matrix(NA, nrow=ndata, ncol=6)

## adjusted rand index between estimated partition and true one
ariTruePartIP <- matrix(NA, nrow=ndata, ncol=6)
ariTruePartCP <- matrix(NA, nrow=ndata, ncol=6)
ariTruePartLSP <- matrix(NA, nrow=ndata, ncol=6)

waicIP <- lmplIP <- matrix(NA, nrow=ndata, ncol=6)
waicCP <- lmplCP <- matrix(NA, nrow=ndata, ncol=6)
waicLSP <- lmplLSP <- matrix(NA, nrow=ndata, ncol=6)

informedPartAlphas <- c(0, 0.25, 0.5, 0.75, 0.9, 0.99)
centeredPartPsis   <- c(0, 10, 20, 50, 80, 100)
lspVariancesTaus <- c(10, 5, 1, 0.05, 1/(N*log(N)), 0.1/(N*log(N)))

colnames(ariInitialPartIP) <- colnames(ariTruePartIP) <- colnames(waicIP) <- colnames(lmplIP) <- as.character(informedPartAlphas)
colnames(ariInitialPartCP) <- colnames(ariTruePartCP)  <- colnames(waicCP) <- colnames(lmplCP) <- as.character(centeredPartPsis)
colnames(ariTruePartLSP) <- colnames(ariTruePartLSP)  <- colnames(waicLSP) <- colnames(lmplLSP) <- as.character(lspVariancesTaus)

####################
## Model settings
####################
## Centered partition 
####################
# Convert data format to use the CP code
## create "intercept" vector
x <- rep(1, N)

## Upper bound for the number of clusters
H_upper <- 50
## number of clusters at initialization
H <- 10

## common priors

## Dirichlet Process concentration parameters
M <- 1
## mean and variance for normal prior on data means
theta <- 0; tau2 <- 10
## format for CP process code
b_prior <-  array(theta, dim = 1)
Q_prior <-  as.matrix(tau2)

########################################################
## Start simulation
########################################################
for(ii in 1:ndata){
  cat("dataset ================================================= ", ii, "\n")
  set.seed(100 + ii);
  cat("seed = ", 100+ii, "\n")

  Ymat <- NULL
  for(t in 1:Tm){
    Ymat <- cbind(Ymat, t(mn[trueP[,t]] + time.dat.gen(N=N, Tm=1, auto_cor=correlation, sd=stdev)))
  if(t == Tm) break
    trueP <- cbind(trueP, c(trueP[,t][-c(1:4)], trueP[,t][c(1:4)]))
  }

  Tmat <- matrix(rep(1:Tm, each=N), nrow=N, byrow=FALSE)

  Yvec <- 1*c(t(Ymat))
  Tvec <- cbind(c(t(Tmat)))
  part.long <- c(t(trueP))

  ## format data for CP 
  Y_list <- as.list(Ymat)
  X_list <- as.list(x)

  Y_list <- lapply(Y_list, function(x) as.matrix(x))
  X_list <- lapply(X_list, function(x) as.matrix(x))

  ##################################################################
  # run our method
  for(j in 1:length(informedPartAlphas)){

    #             # m0, s20, A,            At,   Al,   at, bt, be
    modelPriors <- c(0, 100, 0.5*sd(Yvec), 100,  100,  1,  1,  1)

    cat("informed partition model \n")
    drpm1 <- drpm_fit(y=Ymat, 
              M = 1,    ## dp concentrationa parameter
              initial_partition = initialP, ## initial partition, ## If NULL we have a DP process prior
              ## decide settings for the parameters alpha controlling information from the initial partition
              unit_specific_alpha = FALSE,  ## TRUE - one alpha for each unit
              time_specific_alpha = FALSE,  ## TRUE - one alpha for each time pont
    #                    alphaPriors=rbind(c(1,1)),   # priors for alphas - scalar, or matrix depending on alpha
              starting_alpha = informedPartAlphas[j], 
              alpha_0 = TRUE,               # don't update alpha
              ## model priors - m0, s20, A, 
              modelPriors=modelPriors,
              ## Option for simple model in similations
              simpleModel = 1,  # if 1 then simple model in Sally's simulation
              theta_tau2 =c(0, 10), # this is only used if simpleModel = 1
              ## extra parameters  
    #                   eta1_0 = TRUE, ## not needed
    #                   phi1_0 = TRUE,
              ## MCMC settings
              draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

    cat("centered partition model \n")

    res <- gibbsDensityCP(X_list, Y_list, Z_prior = initialP-1, psi_par = centeredPartPsis[j],
      H, H_upper,alpha = M, b_prior, Q_prior,  niter)
    
    clustSamp = t(res$clustering[, index])
    muSamp = array(0, dim = c(length(Y_list), length(index)))
    for(t in 1:length(index)) {
      muSamp[,t] = res$beta[res$clustering[, index[t]] + 1, ,index[t]]
    }

    cat("LSP partition model \n")
    partitionInit <- sample(1:5, N, replace = T)
    partitionInit <- as.numeric(as.character(factor(partitionInit, labels=order(unique(partitionInit),decreasing=FALSE))))

    resLSP<- LSP_MH(data = as.numeric(Ymat), 
      niter = niterLSP, 
      likelihood_variance = 1 ,
      mu_k = 0, s2_k = 10, 
      tau = lspVariancesTaus[j], 
      initialPart = initialP, 
      partitionInit = partitionInit, 
      tau_proposal = 1/(N*log(N)))

    LSPart <- resLSP$partitionSamples[keepLSP, ]
    muSampLSP <- resLSP$muSamples[keepLSP, ]

    ########## compute output
    cat("compute ARI, LPML, WAIC \n")

    ariMatInitialIP <- ariMatInitialCP <- ariMatInitialLSP <- matrix(NA, nout, Tm)
    ariMatTrueIP <- ariMatTrueCP <- ariMatTrueLSP <-matrix(NA, nout, Tm)

    for(m in 1:nout){
      for(t in 1:Tm) {
        ariMatInitialIP[m,t] <- adjustedRandIndex(drpm1$Si[t,,m], initialP)
        ariMatTrueIP[m,t] <- adjustedRandIndex(drpm1$Si[t,,m], trueP[,t])

        ariMatInitialCP[m,t] <- adjustedRandIndex(clustSamp[m, ], initialP)
        ariMatTrueCP[m,t] <- adjustedRandIndex(clustSamp[m,], trueP[,t])        

        ariMatInitialLSP[m,t] <- adjustedRandIndex(LSPart[m, ], initialP)
        ariMatTrueLSP[m,t] <- adjustedRandIndex(LSPart[m,], trueP[,t])        

      }
    }

    ariInitialPartIP[ii,j] <- mean(apply(ariMatInitialIP,2,mean))
    ariTruePartIP[ii,j] <- mean(apply(ariMatTrueIP,2,mean))

    ariInitialPartCP[ii,j] <- mean(apply(ariMatInitialCP,2,mean))
    ariTruePartCP[ii,j] <- mean(apply(ariMatTrueCP,2,mean))
    
    ariInitialPartLSP[ii,j] <- mean(apply(ariMatInitialLSP,2,mean))
    ariTruePartLSP[ii,j] <- mean(apply(ariMatTrueLSP,2,mean))

    ## GoF for IP prior
    lmplIP[ii, j] <- drpm1$lpml
    waicIP[ii, j] <- drpm1$waic

    ## GoF for CP prior
    xx <- dnorm(unlist(Y_list), muSamp, sd = 1, log = T)
    gofCP <- lpml.robust(xx)
    lmplCP[ii, j] <- gofCP["lpml"]
    waicCP[ii, j] <- gofCP["waic"]

    ## GoF for LSP prior
    xxLSP <- dnorm(unlist(Y_list), muSampLSP, sd = 1, log = T)
    gofLSP <- lpml.robust(xxLSP)
    lmplLSP[ii, j] <- gofLSP["lpml"]
    waicLSP[ii, j] <- gofLSP["waic"]
  

  }
}


saveRDS(list(ariInitialPartIP = ariInitialPartIP, 
             ariTruePartIP = ariTruePartIP,
             ariInitialPartCP = ariInitialPartCP,
             ariTruePartCP = ariTruePartCP,
             ariInitialPartLSP = ariInitialPartLSP,
             ariTruePartLSP = ariTruePartLSP,
             waicIP = waicIP, waicCP = waicCP, 
             lmplIP = lmplIP, lmplCP = lmplCP, 
             lmplLSP = lmplLSP, lmplLSP = lmplLSP,
             scenario = scenario), paste0("output/scenario_",task_id, ".rds"))



