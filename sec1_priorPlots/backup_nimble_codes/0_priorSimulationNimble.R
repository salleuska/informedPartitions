## simulate from prior on temporale partition
## using nimble

library(nimble)
library(mclust) ## adjustedRandIndex
source("RandomPartitionGeneratorNimble.R")

getARI <- function(part1, part2) return(adjustedRandIndex(part1, part2))

## make a function to call R sort function in nimble
getARIR <-  nimbleRcall(prototype   = function(part1 = double(1), part2 = double(1)){}, 
                           Rfun       = "getARI", 
                           returnType = double(0))
   

# computeARI <- nimbleFunction(
# 	run = function( N = double(0),
#   					M = double(0),
#   					rhoVec = double(1),
#   					ntime = double(0), 
#   					FirstPart= double(1)){

# 		out <- array(dim = c(ntime, ntime), type = "double") 
# 		partition <- rtpartitionNimble(N = N, 
# 									   M= M,
# 									   rhoVec = rhoVec,
# 									   ntime = ntime,
# 									   FirstPart = FirstPart)

# 		for(i in 1:ntime){
# 			for(j in 1:ntime){
# 				out[i, j] <- getARIR(partition[, i], partition[, j])

# 			}
# 		}
# 		returnType(double(2))
# 	  	return(out)
# 	}
# )

## Something odd? 
# simulateFunction <- nimbleFunction(
# 	run = function( nRep = double(0),
# 					N = double(0),
#   					M = double(0),
#   					rhoVec = double(1),
#   					ntime = double(0), 
#   					FirstPart= double(1)){

# 		out <- array(dim = c(nRep, ntime, ntime), type = "double") 
# 		partition <- array(dim = c(N, ntime), type = "double") 
# 		for(r in 1:nRep){
# 			partition <<- rtpartitionNimble(N = N, 
# 										   M= M,
# 										   rhoVec = rhoVec,
# 										   ntime = ntime,
# 										   FirstPart = FirstPart)

# 			for(i in 1:ntime){
# 				for(j in 1:ntime){
# 					out[i, j, r] <- getARIR(partition[, i], partition[, j])

# 				}
# 			}

# 		}

# 		returnType(double(3))
# 	  	return(out)
# 	}
# )

# Sim2 
# n of units to partition
N <- 20;
## DP concentration parameter
M <- 0.5;
## time
ntime <- 10;
## Centering partition
FirstPart= seq(1:20)
# FirstPart= rep(c(1,2,3,4), each = 5)
rho <- rep(0.2, N)

# simulate <- simulateFunction(nRep = 10, N = N, M = M, rhoVec = rho, ntime = ntime, FirstPart = FirstPart)
# compileSim <- compileNimble(simulateFunction)

# xx <- compileSim(nRep = 10, N, M, rho, ntime, FirstPart)
simulateP <- rtpartitionNimble(N, M, rho, ntime, FirstPart)
cSimulateP <- compileNimble(simulateP)
# cSimulatePartition(N, M, rho, ntime, FirstPart)
nsim <- 10
out.ours005 <- list()
out.ours025 <- list()
out.ours050 <- list()
out.ours075 <- list()
out.ours095 <- list()
out.ours100 <- list()
# Generate partitions to see differences between our approach and that of Caron
omat7 <- omat8 <- omat9 <- omat10 <- omat11 <- omat12 <- matrix(NA, ntime, ntime)
set.seed(1)
for(ii in 1:nsim){
	if(ii %% 1000 == 0) cat("ii = ", ii, "\n")
	ours005 <- cSimulateP(N=N,M=M,rho=rep(0.05, N),ntime=ntime, FirstPart=FirstPart)
	ours025 <- cSimulateP(N=N,M=M,rho=rep(0.25, N),ntime=ntime, FirstPart=FirstPart)
	ours050 <- cSimulateP(N=N,M=M,rho=rep(0.50, N),ntime=ntime, FirstPart=FirstPart)
	ours075 <- cSimulateP(N=N,M=M,rho=rep(0.75, N),ntime=ntime, FirstPart=FirstPart)
	ours095 <- cSimulateP(N=N,M=M,rho=rep(0.95, N),ntime=ntime, FirstPart=FirstPart)
	ours100 <- cSimulateP(N=N,M=M,rho=rep(1, N),ntime=ntime, FirstPart=FirstPart)

	for(j in 1:ntime){
		for(jj in 1:ntime){

			omat7[j,jj] <- adjustedRandIndex(ours005[,j], ours005[,jj])			
			omat8[j,jj] <- adjustedRandIndex(ours025[,j], ours025[,jj])			
			omat9[j,jj] <- adjustedRandIndex(ours050[,j], ours050[,jj])			
			omat10[j,jj] <- adjustedRandIndex(ours075[,j], ours075[,jj])			
			omat11[j,jj] <- adjustedRandIndex(ours095[,j], ours095[,jj])			
			omat12[j,jj] <- adjustedRandIndex(ours100[,j], ours100[,jj])			
		}
	}
	out.ours005[[ii]] <- omat7
	out.ours025[[ii]] <- omat8
	out.ours050[[ii]] <- omat9
	out.ours075[[ii]] <- omat10
	out.ours095[[ii]] <- omat11
	out.ours100[[ii]] <- omat12
}

