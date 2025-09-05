## need to write an autoregressive process that generate
## alphas = probability of reallocation
####################
## N - number of experimental units
## M - DP scale (dispersion) parameter
## alpha - temporal dependence parameter
## ntime - number of time points
## FirstPart - Allows to provide the first partition if desired
## clusterSpecific - use if partition probability depending on clustering 
####################

logit <- function(x) log(x/(1 - x))
expit <- function(x) exp(x)/(1 + exp(x))

N <- 3;M <- 1;ntime <- 10;FirstPart=c(1,2,1)
# alpha <-  matrix(0.1, ntime, N)
# rtpartition(N, M, alpha, ntime, FirstPart )

## Generate a sequence of alphas (one individual - ntime points)

rep <- 1000
zeta0 <- -3

alpha <- numeric(ntime)
zeta  <- numeric(ntime)

## run simulations multiple times
alphaTest <- matrix(NA, ncol = ntime, nrow = 1000)


for(i in 1:rep){
	kappa  <- runif(1, min = 0, max = 0.5)
	phi     <- rnorm(1, mean = 0, sd = 1)

	epsilon <- rnorm(ntime, mean = 0, sd = kappa)
	zeta[1] <- phi * zeta0 + epsilon[1]
	alpha[1] <- expit(zeta[1])
	for(t in 2:ntime) {
		zeta[t] <-  phi * zeta[t -1] + epsilon[t]
		alpha[t] <- expit(zeta[t])
	}

	alphaTest[i, ] <- alpha

}
hist(alphaTest[, 6], xlim = c(0, 1))

######################
## Test autoregressive partition generator

library(mclust) ## adjustedRandIndex
library(mcclust) ## for variation of information
library(salso) ## for psm - pairwise similarity matrix
###############################
## Source useful functions
###############################
## function that simulate random partitions from the informed partition model
source("utils/randomPartitionGenerator.R")
## F1 distances 
source("utils/F1Distance.R")
###############################
## checks

# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition

FirstPart= rep(c(1,2,3,4), each = 5)
zeta0 <- rep(c(-2.5, -1, 1, 2.5), each = 5)
## stationary process phi in (-1, 1)
## small variance (sd = 0.5)
phi <- -0.5
kappa <- 0.5

nsim <- 5000
clusterSpecific <- FALSE
autoregressive <- TRUE

set.seed(32)

partList <- rtpartition_AR(N = N, M = M, ntime = ntime, 
	FirstPart = FirstPart, 
	kappa = kappa, phi = phi, zeta0 = zeta0 )


library(ggplot2)
library(tidyr)

df <- gather(data.frame(rn = 1:10, t(partList$alphaMat)), key, value, -rn, convert = TRUE) 
df <- df |> filter(key %in% c("X1", "X6", "X11", "X16"))

str(df)	
ggplot(data =df, aes(y = value, group = key, color = key) ) + 
	geom_line(aes(x = rn)) + ylim(c(0,1)) + theme_bw()







