###########
## Functions for F1 measure
f1Measure <- function(cluster1, cluster2, returnLog = FALSE){

 		size1 <- length(cluster1)
 		size2 <- length(cluster2)

 		sizeIntersection <- sum(cluster1 %in% cluster2)

 		f1measureLog <- log(2*sizeIntersection) - log(size1 + size2) 
 		
 		if(returnLog){
 			return(f1measureLog)
 		} else {
 			return(exp(f1measureLog))
 		}
}

# cluster1 <- c(1,2,3)
# cluster2 <- c(1,3,4)
# f1Measure(c(1,2,4), c(1,3,4))

# partition1 <- c(1,1,2,3, 4, 4)
# partition2 <- c(1,1,2,2,3,1)

f1Dist <- function(partition1, partition2, returnLog = FALSE){

	clustSizes1 <- tabulate(partition1)
	clustSizes2 <- tabulate(partition2)
	
	k1 <- length(clustSizes1)
	k2 <- length(clustSizes2)

	intersectionSizes <- table(partition1, partition2)

	allIntNum <- log(2*table(partition1, partition2))
	allIntDen <- log(outer(clustSizes1, clustSizes2, "+"))

	allF1 <- allIntNum - allIntDen

	term1 <- mean(exp(apply(allF1, 1, max)))
	term2 <- mean(exp(apply(allF1, 2, max)))

	out <- 0.5*(term1 + term2)
	if(returnLog){
		log(out)
	} else{
		out
	}
}
