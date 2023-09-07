
######################
## TO DO 
# f1Measure <- nimbleFunction(
#  	run = function(cluster1 = double(1), 
#  		cluster2 = double(1), 
#  		returnLog = TRUE){

#  		size1 <- length(cluster1)
#  		size2 <- length(cluster2)

#  		sizeIntersection <- sum(cluster1 %in% cluster2)

#  		f1measureLog <- log(2*sizeIntersection) - log(size1) - log(size2) 
 		
#  		returnType(double(0))
#  		if(returnLog){
#  			return(f1measureLog)
#  		} else {
#  			return(exp(f1measureLog))
#  		}
# 	}
# )

######################
## LABELS MUST BE SUBSEQUENT
library(nimble)

## There is an error in the relabeling!!!

rtpartitionNimble <- nimbleFunction(
	setup = function(){
	# setup = function(N = double(0),
 #  					M = double(0),
 #  					rhoVec = double(1),
 #  					ntime = double(0), 
 #  					FirstPart= double()){
	# nimCat("check")

	# rhoVec <- numeric(N)
	# if(length(rho) ==  1){
	# 	## SP: if only one number is given, then the pegging probability 
	# 	## is the same for every unit
	# 	rhoVec <- rep(rho, N)
	# } else {
	# 	rhoVec <- rho
	# }

	},
#	run = function(){
 	run = function(N = double(0),
  					M = double(0),
  					rhoVec = double(1),
  					ntime = double(0), 
  					FirstPart= double(1)) {
	# output: n x time matrix of partition
 	# outPart <- array(dim = c(N, ntime), type = "double") 
 	# indMat  <- array(dim = c(N, ntime), type = "double") 


	## output: n x time matrix of partition
 	outPart <- matrix(NA, N, ntime) 
 	indMat  <- matrix(NA, N, ntime)
	
	outPart[, 1] <- FirstPart
	
  clustering <- outPart[, 1]
	indicators <- numeric(N)
	
	## SP: number of clusters for partition at T1
	## current number of clusters
	maxLabel <- max(outPart[, 1])

	clustUniques <- numeric(maxLabel)
	# cluster sizes
  clustCounts <- numeric(maxLabel)
	
  ## count
	aux <-1:maxLabel

  ## Total number of clusters
  K <- 0 
  for(i in seq_along(aux)) { 
    nMembers <- sum(outPart[, 1] == aux[i])
    if(nMembers > 0) {
      clustCounts[aux[i]] <- nMembers   ## cluster size
      clustUniques[K +1] <- aux[i]					## cluster label
      K <- K + 1
    }
	}

	# nimCat("clustCounts  ", clustCounts, "\n")
	# nimCat("clustUniques ", clustUniques, "\n")
 # 	nimCat("K  ", K, "\n")

 # 	nimCat("is there any empty cluster", any(clustCounts == 0))

 	## If there are gaps relabel units
 	##  and update vector of cluster counts
 	if(any(clustCounts == 0)){
	  clustCountsTmp <- numeric(maxLabel)
	  outPartTmp <- outPart[, 1]
  	
  	for(j in 1:K){
  		ind <- outPart[, 1] == clustUniques[j]
  		outPartTmp[ind] <- rep(j, sum(ind)) 
  		clustCountsTmp[j] <- clustCounts[clustUniques[j]]
  	}
  	outPart[, 1] <- outPartTmp
    clustUniques[1:K] <- 1:K
    clustCounts <- clustCountsTmp
  }
 	

  # nimCat("clustCounts after relabeling  ", clustCounts, "\n")
	# nimCat("clustUniques after relabeling", clustUniques, "\n")
 	# nimCat("partition after relabeling  ", outPart[,1], "\n")


  clustering <- outPart[, 1]


	if(ntime > 1) {	 	
		for(t in 2:ntime) { 
			maxLabel <- max(outPart[, t-1])

			subClustCounts <-  numeric(maxLabel)
			subClustUniques <- numeric(maxLabel)
			
			# nimCat("subClustCounts  init ", subClustCounts, "\n")
			# nimCat("subClustUniques init  ", subClustUniques, "\n")

			## sample indicators
			## if 1: resample the observations
			## if 0: peg the observation

			for(i in 1:N){
				indicators[i] <- rbinom(n = 1, size = 1, prob = 1-rhoVec[i])
			}
			
			indSum <- sum(indicators)

			if(indSum == 0){
			## s1: move none (clustering does not change)			
				outPart[, t] <- clustering

			} else if(indSum == N){
			## s2: move all units
			## sample using the DP 
				newClustering <- numeric(N)
				clustCounts <- numeric(N)
				clustUniques <- numeric(N)

				clustCounts[1] <- 1
				clustUniques[1] <- 1
				newClustering[1] <- 1 # cluster label
				K <- 1  # number of clusters
				for(i in 2:N){
					probs <- numeric(K + 1)
					for(k in 1:K){
						probs[k] <- clustCounts[k]/(i-1+M)						
					}	
					probs[K + 1] <- M/(i-1+M)	
					newClustering[i] <- rcat(1, probs)

					## Update set of cluster size and cluster labels
					clustCounts[newClustering[i]] <- subClustCounts[newClustering[i]] + 1
		    	
		    	## if a new cluster has ben sampled, 
		    	## update number of clusters and set of labels
		    	if(newClustering[i] == K + 1){
		    		clustUniques[K + 1] <- K + 1
		    		K <- K + 1
		    	}
		    } 

 			outPart[, t] <- newClustering

			} else {
		## set to 0 the labels of the individuals to re-allocate
			
				## get obs to move and to anchor
				posMove   <- which(indicators == 1)
				posAnchor <- which(indicators == 0)
				

				## set clustering label to 0 for obs. not to move
				clustering[posMove] <- 0
				indMat[, t] <- indicators
				

				subClust <- clustering[posAnchor]
				# # nimCat("subClust ", subClust, "\n")
				
				# # ## Get cluster sizes and number of clusters of
				# # ## the subclustering with allocated units
				aux <- min(subClust):max(subClust)

				# # nimCat("aux subClust round 1 - ", aux, "\n")

		    subK <- 1
		    for(i in seq_along(aux)) { 
		      nMembers <- sum(subClust == aux[i])
		      if(nMembers > 0) {
		        subClustCounts[aux[i]] <- nMembers
		        subClustUniques[subK] <- aux[i]
		        subK <- subK + 1
		      }
		    }
			  subK <- subK - 1

				# # nimCat("subClustUniques round 1 -  ", subClustUniques, "\n")
				# # nimCat("subClustCounts round 1 - ", subClustCounts, "\n")

		    ## relabeling for the subcluster vector  
		    ## if at least one cluster is removed
		    if(subK < K) {
		    	subClustCountsTmp <- numeric(N)

		    	for(j in 1:subK){
		    		subClust[subClust == subClustUniques[j]] <- j 
		    		subClustCountsTmp[j] <- clustCounts[subClustUniques[j]]

		    	}

		    	subClustUniques[1:subK] <- 1:subK
	 	 	  		subClustCounts <- subClustCountsTmp
	 				## update original clustering
	 				clustering[clustering != 0] <- subClust
				}		   
		    # nimCat("subClustUniques end relabeling ", subClustUniques, "\n")
				# nimCat("subClustCounts end relabeling ", subClustCounts, "\n")
	    ##------------------------------------------------------##
	    ## sample units that are not pegged according to the DP
	    ##------------------------------------------------------##
			for(k in seq_along(posMove)){
			#	nimCat("posMove in loop ", posMove[k], "\n")
		    	probs <- numeric(subK + 1)
		    	## compute allocation probability for each subcluster

		    	for(i in 1:subK){
		    		probs[k] <- subClustCounts[k]/(M + sum(subClustCounts)) 
		    	}
		    	probs[subK + 1] <- M/(M + sum(subClustCounts))

					clustering[posMove[k]] <- rcat(1, probs)
					
					## update subclustering and relative counts
					subClust <- clustering[clustering != 0]

					## Update set of cluster size and cluster labels
					subClustCounts[clustering[posMove[k]]] <- subClustCounts[clustering[posMove[k]]] + 1
		    	
		    	## if a new cluster has ben sampled, 
		    	## update number of clusters and set of labels
		    	if(clustering[posMove[k]] == (subK + 1)){
		    		subClustUniques[subK + 1] <- subK + 1
		    		subK <- subK + 1
		    	} 
		    }
				

				## Clustering vector at time T is updated

				aux <-1:max(clustering)

			  ## Total number of clusters
			  K <- 0 
			  for(i in seq_along(aux)) { 
			    nMembers <- sum(outPart[, 1] == aux[i])
			    if(nMembers > 0) {
			      clustCounts[aux[i]] <- nMembers   ## cluster size
			      clustUniques[K + 1] <- aux[i]					## cluster label
			      K <- K + 1
			    }
				}


			 	## If there are gaps relabel units
			 	##  and update vector of cluster counts
			 	 	if(any(clustCounts == 0)){
					  clustCountsTmp <- numeric(maxLabel)
					  clustTmp <- clustering
				  	
				  	for(j in 1:K){
				  		ind <- clustering == clustUniques[j]
				  		clustTmp[ind] <- rep(j, sum(ind)) 
				  		clustCountsTmp[j] <- clustCounts[clustUniques[j]]
				  	}
				  	clustering <- clustTmp
				    clustUniques[1:K] <- 1:K
				    clustCounts <- clustCountsTmp
				  }
				 	

				outPart[, t] <- clustering
			}

		}	
	}
	returnType(double(2))
  return(outPart)
	}
)

N <- 20;M <- 1;rho <- rep(0.5, N);ntime <-20;FirstPart=1:20

# # simulateP <- rtpartitionNimble(N, M, rho, ntime, FirstPart)
# # simulateP$run()
# # cSimulateP <- compileNimble(simulateP, showCompilerOutput = TRUE)

# cSimulateP$run()

## 2
simulateP <- rtpartitionNimble()

# FirstPart=rep(c(1,2,3,4), each = 5)
FirstPart=rep(c(5,2,3,4), each = 5)
ntime <- 10
xx <- simulateP$run(N, M, rho, ntime, FirstPart)


debug(simulateP$run)

# trace(simulateP$run, tracer=quote(if (t==4) browser()))
set.seed(12)
xx <- simulateP$run(N, M, rho, ntime, FirstPart)
table(xx[, 4])

table(xx[,7])

cSimulateP <- compileNimble(simulateP, showCompilerOutput = TRUE)


cSimulateP$run(N, M, rho, ntime, FirstPart)

# undebug(simulateP$run)
# simulateP$run()

# cSimulatePartition <- compileNimble(rtpartitionNimble)

# cSimulatePartition(N, M, rho, ntime, FirstPart)










