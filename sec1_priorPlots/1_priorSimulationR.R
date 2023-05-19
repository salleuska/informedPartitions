library(mclust) ## adjustedRandIndex
library(mcclust) ## for variation of information
library(salso) ## for psm - pairwise similarity matrix

source("RandomPartitionGenerator.R")
source("F1Distance.R")

# cSimulatePartition(N, M, alpha, ntime, FirstPart)
ARIout <- list()
VIout <- list()
F1out <- list()

ARItmp <- VItmp <- F1tmp <- matrix(NA, ntime, ntime)
ARIFirstout<- VIFirstout <- F1Firstout <- matrix(NA, ntime, nsim)
pairWise <- array(NA, dim = c(N, N, ntime, nsim))

set.seed(34)
for(ii in 1:nsim){
	if(ii %% 1000 == 0) cat("ii = ", ii, "\n")
	partList <- rtpartition(N=N,M=M,alpha=alpha,ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)
	
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


settings <- list(alpha = alpha, FirstPart = FirstPart, M = M, N = N, ntime = ntime)
toSave <- list(name = name, 
				settings = settings, 
				pairWise = pairWise, 
				  ARIout = ARIout,
				  VIout = VIout,
				  F1out = F1out,
				  ARIFirstout = ARIFirstout, 
				  VIFirstout = VIFirstout, 
				  F1Firstout = F1Firstout)

saveRDS(toSave, file = paste0("results/samplesDists_firstPart_", name, ".rds"))
