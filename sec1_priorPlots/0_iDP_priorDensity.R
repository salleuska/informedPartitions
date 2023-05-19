#################
## Sally Paganin
## 2023/14/03 
## induced prior on the partition space
#################
library(latex2exp)
library(CPLogit)
source("iCRP_densityFunction.R")
source("makePartitionPlot.R")
alphaVec <- seq(0, 1, length = 11)
####################################
## Example 1 - same pegging probabilities
####################################
## Number of observations
N <- 5;
## DP concentration parameter
M <- 1;

## Centering partition
FirstPart<- c(1,1, 2,2,2) 

## Define scenario for the pegging probabilities
name <- "fixedAlpha"
nVals <- 11
alphaMat <- matrix(NA, nrow = nVals, ncol = length(FirstPart))
alphaMat[, 1:5] <- seq(0, 1, length = nVals)

## generate list of all possible partitions 
partMat <- dist_from(FirstPart - 1, return_partitions = TRUE)$partitions
partMat <- partMat + 1

## compute probabilities 
# tryfoo=function(x){tryCatch(dPDQ(x, part0 = FirstPart, concentration = M, alpha = alphaMat[1, ]),error=function(e){NA})}
# probabilities <- apply(partMat, 1,  function(x) tryfoo(x))

probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$\\alpha_i = \\alpha$")
# plotRes

ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

####################################
## Example 2 - changing probabilities for bigger cluster
####################################
################################
# 2.a cluster 1 probs = 0
FirstPart<- c(1,1,2,2,2) 

## Define scenario for the pegging probabilities
name <- "clust1Alpha0_clust2changing"
nVals <- 11
alphaMat <- matrix(NA, nrow = nVals, ncol = length(FirstPart))
alphaMat[, 1:2] <- 0
alphaMat[, 3:5] <- seq(0, 1, length = nVals)


probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))


plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name, 
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0$")


ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

################################
## 2.a cluster 1 probs = 0.5
alphaMat[, 1:2] <- 0.5

name <- "clust1Alpha05_clust2changing"

probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.5$")

ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

################################
alphaMat[, 1:2] <- 0.6

name <- "clust1Alpha06_clust2changing"
probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.6$")

ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 


################################
FirstPart<- c(1,1,2,2,2) 
alphaMat[, 1:2] <- 0.8

name <- "clust1Alpha08_clust2changing"
probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.8$")

ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

######
FirstPart<- c(1,1,2,2,2) 
alphaMat[, 1:2] <- 0.9

name <- "clust1Alpha09_clust2changing"
probabilities <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.9$")

ggsave(plotRes, file = paste0("iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 


