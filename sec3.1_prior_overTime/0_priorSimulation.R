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

if(!exists("name")) stop("Provide name output file \n")

if(!exists("seedNumber")) {
	seedNumber <- 34
	cat("Using seedNumber ", 34, "\n" )
}
set.seed(seedNumber)

###############################


ARIout <- list()
VIout <- list()
F1out <- list()

ARItmp <- VItmp <- F1tmp <- matrix(NA, ntime, ntime)
ARIFirstout<- VIFirstout <- F1Firstout <- matrix(NA, ntime, nsim)
pairWise <- array(NA, dim = c(N, N, ntime, nsim))

for(ii in 1:nsim){
	if(ii %% 1000 == 0) cat("ii = ", ii, "\n")
	
	if(autoregressive) {
	## Option for autoregressive model
		partList <- rtpartition_AR(N = N, M = M, ntime = ntime, 
			FirstPart = FirstPart, 
			kappa = kappa, phi = phi, zeta0 = zeta0 )
		
	} else {
		partList <- rtpartition(N=N,M=M,alpha=alpha,ntime=ntime, 
			FirstPart=FirstPart, clusterSpecific = clusterSpecific)		
	}
	
	for(j in 1:ntime){
		## pairwise allocation probabilities
		pairWise[,,j, ii] <- psm(partList$ciMat[j, ])
		if(!is.null(FirstPart)){
			ARIFirstout[j,ii] <- adjustedRandIndex(partList$ciMat[j, ], FirstPart)			
			VIFirstout[j,ii] <- vi.dist(FirstPart, partList$ciMat[j, ])/log(N, base = 2)			
			F1Firstout[j,ii] <- f1Dist(FirstPart, partList$ciMat[j, ])						
		}
		for(jj in 1:ntime){
			ARItmp[j,jj] <- adjustedRandIndex(partList$ciMat[j, ], partList$ciMat[jj, ])			
			VItmp[j,jj] <- vi.dist(partList$ciMat[j, ], partList$ciMat[jj, ])/log(N, base = 2)			
			F1tmp[j,jj] <- f1Dist(partList$ciMat[j, ], partList$ciMat[jj, ])			
		}
	}
	ARIout[[ii]] <- ARItmp
	VIout[[ii]] <- VItmp
	F1out[[ii]] <- F1tmp
}


if(autoregressive) {
## Option for autoregressive model
	settings <- list(kappa = kappa, phi = phi, zeta0 = zeta0, 
		FirstPart = FirstPart, M = M, N = N, ntime = ntime)
	
} else {
	settings <- list(alpha = alpha, FirstPart = FirstPart, M = M, N = N, ntime = ntime)
}


toSave <- list(name = name, 
				settings = settings, 
				pairWise = pairWise, 
				  ARIout = ARIout,
				  VIout = VIout,
				  F1out = F1out,
				  ARIFirstout = ARIFirstout, 
				  VIFirstout = VIFirstout, 
				  F1Firstout = F1Firstout, 
				  seedNumber = seedNumber)

saveRDS(toSave, file = paste0("results/samples_", name, ".rds"))
