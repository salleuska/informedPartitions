############################
## Function for relabeling
relabel <- function(partition){
	out <- as.numeric(factor(partition, levels = unique(partition)))
	return(out)
}

## Dirchlet process EPPF 
dpEPPF <- function(partition, concentration, returnLog = FALSE){
  require(betafunctions)
  nObs = length(partition)
  nBlocks = length(unique(partition))
  blockSizes = table(partition)
  
  logProb = nBlocks*log(concentration) + sum(lgamma(blockSizes)) - log(afac(concentration, nObs))
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}


############################

# gamma <- c(0,0,1)
# part0 <- c(0,0,1)
# partition <- c(0,1,2)
# concentration <- 1

## conditional partition distribution for the 
## informed Dirchlet process (IDP)

dIDPConditional <- function(partition, part0, concentration, gamma) {
	require(CPLogit)

	## This computes the probability of a partition conditional to part0 and gamma

	partition = relabel(partition)
	part0 = relabel(part0)

	m <- length(partition)
	## Generate all possible partitions
	setPartSpace <- dist_from(rep(0, m), return_partitions = TRUE)$partitions
	nPart <- dim(setPartSpace)[1]

	## compute partition probability for all compatible ones
	## 1. check compatibility; if partition is not compatible then the prob = 0

	## derive reduced partitions (relabeled)
	reducedPart0 <- relabel(part0[as.logical(gamma)])
	reducedPartition <- relabel(partition[as.logical(gamma)])

	compatibility <- all(reducedPart0 == reducedPartition)

	if(compatibility) {
		## find all partitions in the space compatible with part0 and gamma
		## 1. calculate the reduced partition space
		reducedPartSpace    <- setPartSpace[, as.logical(gamma), drop=FALSE]
		reducedPartSpace[]  <- t(apply(reducedPartSpace, 1, relabel))
		## 2. compatibility checks with reduced part0 and reducedPartSpace
		compatiblePartRows <- apply(reducedPartSpace, 1, function(x) all(reducedPart0 == x))
		compatibleSetPart  <- setPartSpace[compatiblePartRows,,drop = FALSE]
		
		logProb <- dpEPPF(partition, concentration = concentration, returnLog = TRUE) 
		logNormConst <- log(sum(apply(compatibleSetPart, 1, function(x) dpEPPF(x, concentration = concentration,  returnLog = FALSE))))

		prob <- exp(logProb - logNormConst)

	} else {
		prob <- 0
	}

	return(prob)


}

# gamma <- c(1,1,1)
# part0 <- c(0,1,2)
# partition <- c(0,1,2)
# concentration <- 1
# alpha <- 0.1

# dIDPConditional(partition, part0 , concentration, gamma)

## density function for the informed Dirichlet Process

dIDP <- function(partition, part0, concentration, alpha){
	partition <- relabel(partition)
	part0 <- relabel(part0)

	m <- length(partition)
	if(length(alpha) == 1) 	alphaVec <- rep(alpha, m) else if (length(alpha) == m) alphaVec <- alpha else stop("Error: alpha needs to be a vector of size 1 or ", m, "\n")

	gammaSpace <- as.matrix(expand.grid(replicate(m, 0:1, simplify = FALSE)))
	nGamma <- dim(gammaSpace)[1] 

	gammaProbs <- numeric(m)
	gammaProbs <- t(sapply(1:NROW(gammaSpace), \(i) alphaVec^gammaSpace[i,] * (1-alphaVec)^(1 - gammaSpace[i,])))

	probsVecTmp <- numeric(m)
	for(i in 1:nGamma){
		## check if partition is compatible with part0 according to the gamma Vector
		reducedPart0 <- relabel(part0[as.logical(gammaSpace[i, ])])
		reducedPartition <- relabel(partition[as.logical(gammaSpace[i, ])])

		compatibility <- all(reducedPart0 == reducedPartition)		
		if(compatibility){ 
			probsVecTmp[i] <- dIDPConditional(partition, part0, concentration, gamma = gammaSpace[i, ])
		} else probsVecTmp[i] <- 0
	}

	## probability
	logComponents <- sapply(1:nGamma, function(x) sum(log(gammaProbs[x, ])) + log(probsVecTmp[x])) 
	
	return(sum(exp(logComponents)))
}

# dIDP(partition, part0, concentration =  1, alpha = c(0.2, 0.1, 0.1))
## The next two should be equivalent
# dPDQ(partition, part0, concentration =  1, alpha = 0)
# dpEPPF(partition,concentration =  1)





