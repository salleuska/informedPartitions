#################
## 2022/12/01 
## induced prior on the partition space
#################
library(mclust) ## adjustedRandIndex
library(mcclust) ## for variation of information
library(salso) ## for psm - pairwise similarity matrix
library(latex2exp)
library(CPLogit)
source("RandomPartitionGenerator.R")
source("makePartitionPlot.R")
#################
N <- 5;
## DP concentration parameter
M <- 1;
## time
ntime <- 2;
## Centering partition
# FirstPart= seq(1:20)
name <- "fixedAlpha"
FirstPart<- c(1,1, 2,2,2) 

alphaVec <- seq(0, 1, length = 11)

nsim <- 5000
clusterSpecific <- FALSE

partMat <- dist_from(FirstPart - 1, return_partitions = TRUE)$partitions
partMat <- partMat + 1
partVec <- apply(partMat, 1, function(x)  paste0(x, collapse = ""))

# partVec <- c("12313","11111","12121","11212","12113","12222",
# 	"11211","12132","12311","12112","12111","12212",
# 	"12342","11213","12343","11121","11221","11122",
# 	"12134","12333","12344","12345","12221","11112",
# 	"12223","12131","12321","12322","12232","11231",
# 	"12211","12123","12324","12233","12213","12332",
# 	"12331","12341","11123","11222","11234","12323",
# 	"11232","12314","11233","12234","11223","12122",
# 	"12334","12312","12133","12231")

# library(stringr)
# partMat <- str_split_fixed(partVec, pattern = "", n = nchar(partVec))
# partMat <- apply(partMat, 2, as.numeric)

probs <- matrix(NA, nrow = bell(5), ncol = length(alphaVec))
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaVec[i],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name, 
			xLabel = "$\\alpha_i = \\alpha$")

ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 



################################
name <- "clust1Alpha0_clust2changing"
FirstPart<- c(1,1,2,2,2) 

nVals <- 11

alphaMat <- matrix(NA, nrow = nVals, ncol = length(FirstPart))
alphaMat[, 1:2] <- 0
alphaMat[, 3:5] <- seq(0, 1, length = nVals)

nsim <- 5000
clusterSpecific <- FALSE


probs <- matrix(NA, nrow = bell(5), ncol = nVals)
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaMat[i, ],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name, 
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0$")


ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

################################
alphaMat[, 1:2] <- 0.5

name <- "clust1Alpha05_clust2changing"

probs <- matrix(NA, nrow = bell(5), ncol = nVals)
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaMat[i, ],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.5$")

ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 
################################
FirstPart<- c(1,1,2,2,2) 
alphaMat[, 1:2] <- 0.6

name <- "clust1Alpha06_clust2changing"

probs <- matrix(NA, nrow = bell(5), ncol = nVals)
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaMat[i, ],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.6$")

ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 


################################
FirstPart<- c(1,1,2,2,2) 
alphaMat[, 1:2] <- 0.8

name <- "clust1Alpha08_clust2changing"

probs <- matrix(NA, nrow = bell(5), ncol = nVals)
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaMat[i, ],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.8$")

ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 

######
FirstPart<- c(1,1,2,2,2) 
alphaMat[, 1:2] <- 0.9

name <- "clust1Alpha09_clust2changing"

probs <- matrix(NA, nrow = bell(5), ncol = nVals)
for(i in 1:length(alphaVec)){
	set.seed(134 + i*2)
	out <- replicate(nsim, rtpartition(N=N,M=M,alpha=alphaMat[i, ],ntime=ntime, FirstPart=FirstPart, clusterSpecific = clusterSpecific)$ciMat)

	parts <- apply(out[2, ,], 2, function(x) paste0(x, sep = "", collapse = ""))


	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}

saveRDS(list(probs = probs, c0 = FirstPart, partVec = partVec), file = paste0("partitionsDensity_", name, ".rds"))

plotRes <- makePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = name,
			xLabel = "$(\\alpha_{3}, \\alpha_{4}, \\alpha_{5}) = \\alpha$ for fixed $\\alpha_{1} = \\alpha_{2} = 0.9$")

ggsave(plotRes, file = paste0("SCDRPM_", name, ".pdf"), 
		device = "pdf", height = 20, width = 25, unit = "cm", dpi = 300) 




