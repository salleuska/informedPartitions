library(partitions)
########################
## allParts <- setparts(5)
allParts <- CPLogit::dist_from(c(0,0,0,0,0), 
	return_partitions = TRUE)$partitions + 1

allParts

########################
## Sequential allocation for DP
## conditional to one permutation
## alpha : concentration parameter
## n: number of elements

# dSequentialDP : density
# rSequentialDP : sampling

## convention that first element is allocated to first st
## labels are ordered

dSequentialDP <- function(part, alpha, returnLog = FALSE){
	## conditional to a fixed permutation
	n <- length(part)
	logProb <- numeric(n)

	logProb[1] <- 0
	nClust <- 1

	for(j in 2:n) {
		nClustNew <- length(unique(part[1:j])) 
		if(nClustNew == nClust){
			## element j is in one of the occupied clusters
			nClustj <- sum(part[1:j] == part[j]) - 1
			logProb[j] <- log(nClustj) - log(alpha + j - 1)


		} else {
			## element j is in a new cluster
			logProb[j] <- log(alpha) - log(alpha + j - 1)
		}
		nClust <- nClustNew
	}

	if(returnLog){
		sum(logProb)
	} else{
		exp(sum(logProb))
	}
}

## As is, this is already normalized
# sum(apply(allParts, 1, function(c) dSequentialDP(c, 1)))

# probs <- apply(allParts, 1, function(x) dSequentialDP(x, 1))
# probs <- apply(allParts, 1, function(x) dSequentialDP(x, 1, returnLog = T))


dSeqDPmarginal <- function(part, alpha){
	require(gtools)
	n <- length(part)
	perm <- permutations(n = n, r = n, v = 1:n)     

	nPermutations <- dim(perm)[1]

	prob <- sapply(1:nPermutations, function(x) 
		dSequentialDP(part[perm[x, ]], alpha = alpha))
	prob
}


probsMarginalTmp <- t(apply(allParts, 1, function(x) dSeqDPmarginal(x, 1)))
probsMarginal <- apply(probsMarginalTmp, 1, sum)
probsMarginal <- probsMarginal/sum(probsMarginal)
sum(probsMarginal)

####################
## output of dSequetial DP should match EPPF 
dpEPPF <- function(partition, concentration, returnLog = FALSE){
  nObs = length(partition)
  nBlocks = length(unique(partition))
  blockSizes = table(partition)
  
  logProb = nBlocks*log(concentration) + sum(lgamma(blockSizes))
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}

probs2 <- apply(allParts, 1, function(x) dpEPPF(x, 1))
probs2 <- probs2/sum(probs2)

plot(probsMarginal, probs2)
abline(0, 1)


##### DD package - shallot ##################
# library(shallot)
# p <- ewens(mass(1), n.items = 3)
# dpEPPFdahl <- partition.pmf(ewens(mass(1), n.items = 3))

# probsDahl <- apply(allParts, 1, function(x) dpEPPFdahl(x, log = FALSE))
# probsDahl
#############################################

### Discrepancy between subsets - F1 measure
## Gentz et al
# part1 <- c(1,1,2,3,2)
# part2 <- allParts[3,]
# ## 1 1 1 2 3

# unit <- 3

# sdiscF1(part1, part2, 1)
# sdiscF1(part1, part2, 2)
# sdiscF1(part1, part2, 3)
# sdiscF1(part1, part2, 5)


sdiscF1 <- function(part1, part2, unit, returnLog = FALSE) {
	
	sub1Label <- part1[unit]	
	sub2Label <- part2[unit]	

	sub1Size <- sum(part1[1:unit] == part1[unit])
	sub2Size <- sum(part2[1:unit] == part2[unit])

	intersectionSize <- length(intersect(which(part1[1:unit] == part1[unit]), 
														which(part2[1:unit] == part2[unit])))

	logDisc <- log(sub1Size + sub2Size - 2*intersectionSize) - 
			log(sub1Size + sub2Size)

	if(returnLog){
		logDisc
	} else{
		exp(logDisc)
	}
}
#######################
part <- allParts[11,]
alpha <- 1

basePartition <- c(1,2,1,2,3)

dSeqDPCentered <- function(part, alpha, 
						basePartition,
						discrepancy = "sdiscF1", 
						psi = 1,
						returnLog = FALSE){
	## conditional to a fixed permutation
	n <- length(part)
	logProb <- numeric(n)
	contributionCentering <- numeric(n)

	logProb[1] <- 0
	nClust <- 1

	for(j in 2:n) {
		nClustNew <- length(unique(part[1:j])) 
		if(nClustNew == nClust){
			## element j is in one of the occupied clusters
			nClustj <- sum(part[1:j] == part[j]) - 1
			contributionDP <- log(nClustj) - log(alpha + j - 1)

		} else {
			## element j is in a new cluster
			contributionDP <- log(alpha) - log(alpha + j - 1)
		}

		contributionCentering[j] <- do.call(discrepancy, list(part1  =part, 
																		 part2 = basePartition, 
																		 unit = j, 
																		 returnLog = FALSE))
		
		logProb[j] <- contributionDP - psi * contributionCentering[j]		

		nClust <- nClustNew
	}

	normConst <- log(cumsum(exp(-psi*contributionCentering)))

	if(returnLog){
		sum(logProb - normConst)
	} else{
		exp(sum(logProb - normConst))
	}
}

dSeqDPCenteredMarginal <- function(part, alpha, basePartition, discrepancy, psi){
	require(gtools)
	n <- length(part)
	perm <- permutations(n = n, r = n, v = 1:n)     

	nPermutations <- dim(perm)[1]

	prob <- sapply(1:nPermutations, function(x) 
		dSeqDPCentered(part[perm[x, ]], 
			alpha = alpha, 
			basePartition = basePartition, 
			psi = psi, 
			discrepancy = discrepancy))
	prob
}w

basePartition <- c(1,2,3,4,4)
probs <- t(apply(allParts, 1, function(x) 
				dSeqDPCenteredMarginal(x, 
											  alpha = 1, 
											  psi = 2,
											  basePartition =  basePartition, 
											  discrepancy = "sdiscF1")))
probsMarginal <- apply(probs, 1, sum)
probsMarginal <- probsMarginal/sum(probsMarginal)
probsMarginal 
