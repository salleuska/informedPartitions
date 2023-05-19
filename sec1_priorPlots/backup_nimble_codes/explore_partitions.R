# library(partitions)
# ########################
# ## allParts <- setparts(5)
# allParts <- CPLogit::dist_from(c(0,0,0,0,0), 
# 	return_partitions = TRUE)$partitions + 1
# allParts
library(mcclust) ## for variation of information
library(mclust) ## adjustedRandIndex
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

# 1 - f1Dist(partition1, partition2)
# vi.dist(partition1, partition2)/log(length(partition2), base = 2)



#########################
source("RandomPartitionGenerator_sal.R")
#########################
## n of units to partition
# N <- 5;
# ## DP concentration parameter
# M <- 1;
# ## time
# ntime <- 10;
# ## Centering partition
# FirstPart=c(1,1,1, 2,2)
##########
# Sim2 
# n of units to partition
N <- 20;
## DP concentration parameter
M <- 1;
## time
ntime <- 10;
## Centering partition
FirstPart= rep(c(1,2,3,4), each = 5)
####################################

# ## Function for when I want to simulate alpha from a beta prior
# simPart <- function(N, M, ntime, FirstPart){
# 	rho <- matrix(0, nrow = ntime, ncol = N)
# 	clustSizes <- tabulate(FirstPart)

# 	## simulate pegging probabilities
# 	for(i in 1:length(clustSizes)){
# 		rho[, 1:clustSizes[i]] <- rbeta(clustSizes[i]*ntime, clustSizes[i], 1)
# 	}
# 	clustMat <- rtpartition(N, M, rho, ntime, FirstPart)$ciMat
# 	clustMat 
# }
# set.seed(1234)

# matPart <- replicate(10000, simPart(N, M, ntime, FirstPart))

rho <- matrix(rep(c(0.2,0.4,0.6,0.8), each = 5*10), nrow = ntime, ncol = N)
set.seed(4354)
matPart <- replicate(10000, rtpartition(N, M, rho, ntime, FirstPart)$ciMat)

distVi <- apply(matPart, c(1, 3), function(x) vi.dist(FirstPart, x))
distVi <- distVi/log(N, base = 2)

distF1 <- apply(matPart, c(1, 3), function(x) 1- f1Dist(FirstPart, x))
# distF1 <- 1 - distF1
settings <- list(N = N, M = M, ntime = ntime, FirstPart = FirstPart, rho =rho,
				 matPart = matPart, distVi = distVi, distF1 = distF1)	
save(settings, file = "simN20_rhocluster.rds")
###########
library(ggplot2)
library(ggridges)
theme_set(theme_minimal())
library(dplyr)
library(tibble)
library(ClusterR)

distDF <- reshape2::melt(distF1)

distDF = as_tibble(distDF) |>
	mutate(time = factor(Var1, levels =  1:length(unique(distDF$Var1))))

str(distDF)

ggplot(distDF) +
    geom_density_ridges(aes(x = value, y = time))
ggplot(distDF) +
    geom_density_ridges(aes(x = value, y = time),
    stat = "binline", bins = 100, scale = 0.95,
    draw_baseline = FALSE)



############
## function for time-lagged partiton/clsuter distances

randMat <- matrix(NA, ntime, ntime)
viMat <- matrix(NA, ntime, ntime)
f1Mat <- matrix(NA, ntime, ntime)
## matrix for initiaal partition subsets? 

smatPart[1,2,]


matPart[,,1]

# 1 - f1Dist(matPart[1,,4], matPart[5,,4])

# for(i in 1:ntime){
# 	cat("i = ", i, "\n")
# 	for(j in i:ntime){
# 		f1Mat[i, j] <- mean(apply(matPart, 3, function(x) 1- f1Dist(x[i,], x[j,])))
# 	}
# }

str(f1Mat)





