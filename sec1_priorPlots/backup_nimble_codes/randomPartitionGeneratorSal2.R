## Functions needed to generate data

# November 2 2018
# N - number of experimental units
# M - DP scale (dispesrsion) parameter
# rho - temporal dependence parameter
# ntime - number of time points
# tau - sd associated with atom generation
# sig - sd associated with data generation
# Caron - logical indicating if temporal dependence follows Caron
# FirstPart - Allows me to provide the first partition if so desired
# Space - logical indicating of sPPM should be used to generate partition or not

# N <- 3;M <- 1;rho <- 0.5;ntime <- 2;Caron=FALSE;FirstPart=NULL

# FirstPart <- c(1,2,3)

rtpartition <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL){
	library(plyr)
	out <- NULL

	clustering <- FirstPart

	if(is.null(FirstPart)){
		# Create cluster at time period t = 1
		clustering <- 1 # cluster label
		clustSizes <- 1 # vector of cluster sizes
		K <- 1 # number of clusters

		for(k in 2:N){

			p <- c(clustSizes/(k-1+M), M/(k-1+M))
			clustering <- c(clustering, sample(1:(K+1), 1, prob=p))
			clustSizes <- table(clustering)
			K <- length(unique(clustering))
		}

	}

#	print(clustering)
#	print(mustar)
#	plot(1:N, clustering, col=clustering, pch=19, cex=1.25)

	# Now to generate 10-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	clusteringMat <- matrix(NA, nrow=ntime, ncol=N)
	clusteringMat[1,]  <- clustering
	K <-  length(unique(clustering))
	indicatorMat <- matrix(0,nrow=ntime, ncol=N)

#	print(clusteringMat)
	if(ntime > 1){	
		for(t in 2:ntime){
#				cat("t = ", t, "\n")
			nObsToSample <- rbinom(1, N, 1-rho)
			obsToSample <- sample(1:length(clustering), nObsToSample)
			if(nObsToSample > 0){
				clustering[obsToSample]<- 0
				indicatorMat[t,obsToSample] <- 1

			}	

#				print(clustering)

			clustSizes <- tabulate(clustering[clustering!=0]);
			if((K - length(clustSizes)) > 0){
				clustSizes <- c(clustSizes, rep(0, K-length(clustSizes)))
			}

			K <- length(unique(clustering[clustering!=0]))


#				print(K); print(length(clustSizes))

#				print(clustSizes)

			if(r < N){
			
				clustering[clustering!=0]<-as.numeric(as.character(factor(clustering[clustering!=0], labels=order(unique(clustering[clustering!=0]),decreasing=FALSE))))

			}
#				print(clustering)



			for(k in dnk){
#					cat("k = ", k, "\n")
#					print(clustSizes)
				p <- 1
				if(K>0) p <- c(clustSizes[clustSizes!=0]/(sum(clustSizes[clustSizes!=0])+M), M/(sum(clustSizes[clustSizes!=0])+M))
#					print(p)
#					cat("K + 1 = ", K+1, "\n")

				clustering[k] <- sample(1:(K+1), 1, prob=p)
				clustSizes <- table(clustering[clustering!=0])
				K <- length(unique(clustering[clustering!=0]))
#					print(clustering)
#					print(table(clustering))
			}
#				print(clustering)
			clusteringtmp <- factor(clustering,  labels=order(unique(clustering), decreasing=FALSE))

#				print(clusteringtmp)
			clusteringMat[t,] <- as.numeric(as.character(clusteringtmp))
#				print(clusteringMat)

		}
	
	}
	out$clusteringMat <- clusteringMat
	out$indicatorMat <- indicatorMat
	out
}
