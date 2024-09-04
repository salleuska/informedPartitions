####################
## Function to simulate partitions with dependence 
## sover time from the prior distribution
####################
## N - number of experimental units
## M - DP scale (dispersion) parameter
## alpha - temporal dependence parameter
## ntime - number of time points
## FirstPart - Allows to provide the first partition if desired
## clusterSpecific - use if partition probability depending on clustering 

# N <- 3;M <- 1;ntime <- 10;FirstPart=c(1,2,1)
# alpha <-  matrix(0.1, ntime, N)
# rtpartition(N, M, alpha, ntime, FirstPart )

isScalar <- function(x) is.atomic(x) && length(x) == 1L
isVector <- function(x) is.atomic(x) && is.vector(x) 
	
rtpartition <- function(N, M ,alpha,ntime, FirstPart=NULL, 
	clusterSpecific = FALSE){
	library(plyr)
	out <- NULL

	if(isScalar(alpha)){
		message("using same alpha for each observation over time")
		alphaMat <- matrix(alpha, nrow = N, ncol = ntime)
	} else if(isVector(alpha)) {
		message("using constant alpha across time")
		alphaMat <- matrix(rep(alpha, ntime), nrow = N, ncol = ntime)
	} else if(is.matrix(alpha)){
		if(clusterSpecific){
			## if not using the cluster dependent probability
			## then behave as the DP 
			alphaMat <- matrix(0, nrow = N, ncol = ntime)
			alphaMat[1:dim(alpha)[1], 1:dim(alpha)[2]] <- alpha
		} else{
			alphaMat <- alpha
		}
	}

	## initialize clustering with a given partition
	clustering <- FirstPart

	## SP: if the initial partition is not given
	## then simulate the partition from the DP
	if(is.null(FirstPart)){
		# sample  cluster at time period t = 1
		clustering <- 1 # cluster label
		clustSizes <- 1 # vector of cluster sizes
		K <- 1  # number of clusters

		for(k in 2:N){

			p <- c(clustSizes/(k-1+M), M/(k-1+M))
			clustering <- c(clustering, sample(1:(K+1), 1, prob=p))
			clustSizes <- table(clustering)
			K <- length(unique(clustering))
		}

	}
   
#   print(clustering)
#	print(mustar)
#	plot(1:N, clustering, col=clustering, pch=19, cex=1.25)

	# Now to generate ntime-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	
	clusteringMat <- matrix(NA, nrow=ntime, ncol=N) ## matrix of partition sequences
	clusteringMat[1,]  <- clustering

	## SP: number of clusters for partition at T1
	K <-  length(unique(clustering))
	gMat <- matrix(0,nrow=ntime, ncol=N) ## matrix of pegging vector indicators

#	print(clusteringMat)
	if(ntime > 1){
		
		for(t in 2:ntime){
#		cat("t = ", t, "\n")
			if(clusterSpecific){
				indicatorVec <- rbinom(N, 1, 1 - alphaMat[clustering, t])
			}	else {
				indicatorVec <- rbinom(N, 1, 1 - alphaMat[, t])
			}		
			r <- sum(indicatorVec)
			obsToSample <- which(indicatorVec == 1)
			if(r > 0){
				clustering[obsToSample] <- 0
				gMat[t, ] <- indicatorVec

			}	


#			print(clustering)

			## SP: get the cluster sizes for labels from 1 to K
			clustSizes <- tabulate(clustering[clustering!=0])

			if((K - length(clustSizes)) > 0){
				clustSizes <- c(clustSizes, rep(0, K-length(clustSizes)))
			}
		
			K <- length(unique(clustering[clustering!=0]))


#				print(K); print(length(clustSizes))

#				print(clustSizes)

			if(r < N){
				## relabeling of observations
				clustering[clustering!=0] <- as.numeric(as.character(factor(clustering[clustering!=0], labels=order(unique(clustering[clustering!=0]),decreasing=FALSE))))

			}
#				print(clustering)



			for(k in obsToSample){
#				cat("k = ", k, "\n")
#				print(clustSizes)
				allocationProb <- 1
				if(K>0) allocationProb <- c(clustSizes[clustSizes!=0]/(sum(clustSizes[clustSizes!=0])+M), M/(sum(clustSizes[clustSizes!=0])+M))
				
				## SP: p is the vector of allocation probabilities to the 
				## already occupied clusters + new cluster
				# print(allocationProb)
				# cat("K + 1 = ", K+1, "\n")

				## update partition
				clustering[k] <- sample(1:(K+1), 1, prob=allocationProb)
				## update cluster sizes and number of clusters 
				clustSizes <- table(clustering[clustering!=0])
				K <- length(unique(clustering[clustering!=0]))
#					print(clustering)
#					print(table(clustering))
			}
#			print(clustering)
			## relabeling of observations (no gaps)
			clusteringtmp <- factor(clustering,  labels=order(unique(clustering), decreasing=FALSE))

#			print(clusteringtmp)
			clusteringMat[t,] <- as.numeric(as.character(clusteringtmp))
#			print(clusteringMat)

		}
		
	}
	
	out$ciMat <- clusteringMat
	out$gMat <- gMat
	out
}

####################################################
## Generate partitions usign latent autoregressive 
## process for the alphas 
####################################################
# N <- 3;M <- 1;ntime <- 10;FirstPart=c(1,2,1)
# alpha <-  matrix(0.1, N, ntime)
# rtpartition(N, M, alpha, ntime, FirstPart )

# phi_sd = 1
# kappaMin = 0.01
# kappaMax = 0.5
# zet0 <- 3

# rtpartition_AR(N = N, M = M, ntime = ntime, FirstPart = FirstPart, 
# 	kappa = 0.5, phi = 1, zeta0 = 3)

rtpartition_AR <- function(N, M, ntime, 
	## first partition
	FirstPart=NULL, 
	## first coeff
	zeta0,
	kappa = 0.5, 
	phi = 1,  
	## priors hyperameters	
	clusterSpecific = FALSE,
	hyperpriors = FALSE, 
	hyperparameters = NULL){

	library(plyr)
	out <- NULL

    ## hyperprior parameters in case I want to add a layer of simulation
	if(isScalar(zeta0)){
		message("using same zeta0 for each observation")
		zeta0 <- rep(zeta0, N)
	}
	## initialize matrices
	alphaMat <- matrix(NA, nrow = N, ncol = ntime)
	zetaMat <- matrix(NA, nrow = N, ncol = ntime)
	epsilon <- matrix(NA, nrow = N, ncol = ntime)

	if(hyperpriors) {
		hyperparameters$phi_sd <-  1
		hyperparameters$kappaMin <-  0.01
		hyperparameters$kappaMax <-  1
		kappas <- runif(N, min = kappaMin, max = kappaMax)
		phi    <- rnorm(N, mean = 0, sd = phi_sd)
	} else {
		kappas <- rep(kappa, N) 
		phi    <- rep(phi, N)

	}

	zetaMat[, 1] <- phi * zeta0 + rnorm(N, mean = 0, sd = kappas)
	alphaMat[, 1] <- expit(zetaMat[,1])
	for(t in 2:ntime) {
		zetaMat[, t] <-  phi * zetaMat[, t -1] +  rnorm(N, mean = 0, sd = kappas)
		alphaMat[, t] <- expit(zetaMat[, t])
	}

	## initialize clustering with a given partition
	clustering <- FirstPart

	## SP: if the initial partition is not given
	## then simulate the partition from the DP
	if(is.null(FirstPart)){
		# sample  cluster at time period t = 1
		clustering <- 1 # cluster label
		clustSizes <- 1 # vector of cluster sizes
		K <- 1  # number of clusters

		for(k in 2:N){

			p <- c(clustSizes/(k-1+M), M/(k-1+M))
			clustering <- c(clustering, sample(1:(K+1), 1, prob=p))
			clustSizes <- table(clustering)
			K <- length(unique(clustering))
		}

	}
   
#   print(clustering)
#	print(mustar)
#	plot(1:N, clustering, col=clustering, pch=19, cex=1.25)

	# Now to generate ntime-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	
	clusteringMat <- matrix(NA, nrow=ntime, ncol=N) ## matrix of partition sequences
	clusteringMat[1,]  <- clustering

	## SP: number of clusters for partition at T1
	K <-  length(unique(clustering))
	gMat <- matrix(0,nrow=ntime, ncol=N) ## matrix of pegging vector indicators

#	print(clusteringMat)
	if(ntime > 1){
		
		for(t in 2:ntime){
#		cat("t = ", t, "\n")
			if(clusterSpecific){
				indicatorVec <- rbinom(N, 1, 1 - alphaMat[clustering, t])
			}	else {
				indicatorVec <- rbinom(N, 1, 1 - alphaMat[, t])
			}		
			r <- sum(indicatorVec)
			obsToSample <- which(indicatorVec == 1)
			if(r > 0){
				clustering[obsToSample] <- 0
				gMat[t, ] <- indicatorVec

			}	


#			print(clustering)

			## SP: get the cluster sizes for labels from 1 to K
			clustSizes <- tabulate(clustering[clustering!=0])

			if((K - length(clustSizes)) > 0){
				clustSizes <- c(clustSizes, rep(0, K-length(clustSizes)))
			}
		
			K <- length(unique(clustering[clustering!=0]))


#				print(K); print(length(clustSizes))

#				print(clustSizes)

			if(r < N){
				## relabeling of observations
				clustering[clustering!=0] <- as.numeric(as.character(factor(clustering[clustering!=0], labels=order(unique(clustering[clustering!=0]),decreasing=FALSE))))

			}
#				print(clustering)



			for(k in obsToSample){
#				cat("k = ", k, "\n")
#				print(clustSizes)
				allocationProb <- 1
				if(K>0) allocationProb <- c(clustSizes[clustSizes!=0]/(sum(clustSizes[clustSizes!=0])+M), M/(sum(clustSizes[clustSizes!=0])+M))
				
				## SP: p is the vector of allocation probabilities to the 
				## already occupied clusters + new cluster
				# print(allocationProb)
				# cat("K + 1 = ", K+1, "\n")

				## update partition
				clustering[k] <- sample(1:(K+1), 1, prob=allocationProb)
				## update cluster sizes and number of clusters 
				clustSizes <- table(clustering[clustering!=0])
				K <- length(unique(clustering[clustering!=0]))
#					print(clustering)
#					print(table(clustering))
			}
#			print(clustering)
			## relabeling of observations (no gaps)
			clusteringtmp <- factor(clustering,  labels=order(unique(clustering), decreasing=FALSE))

#			print(clusteringtmp)
			clusteringMat[t,] <- as.numeric(as.character(clusteringtmp))
#			print(clusteringMat)

		}
		
	}
	
	out$ciMat <- clusteringMat
	out$gMat <- gMat
	out
}

