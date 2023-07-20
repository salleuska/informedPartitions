#################
## Sally Paganin
## 2023/14/03 
#################
## This script produces a density plot over the set partition space of 5 elements
## for the LSP model (smith and allenby 2021)
#################
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(reshape2)

library(CPLogit)
require(grid)
require(gridExtra)

## original
# sourceCpp("dmrpfunctions.cpp")
sourceCpp("LSP/dmrpfunctions.cpp")
source("makeBasePartitionPlot.R")
####################################
## values of the scale parameter tau
tauVec <- seq(3, 0, length =11)
####################################
## Number of observations
N <- 5;
## DP concentration parameter
M <- 1;

## Centering partition
FirstPart<- c(1,1, 2,2,2) 


partMat <- dist_from(FirstPart - 1, return_partitions = TRUE)$partitions
partMat <- partMat + 1
partVec <- apply(partMat, 1, function(x)  paste0(x, collapse = ""))

set.seed(1)
nsim <- 30000

rho <- FirstPart
n <- length(FirstPart)

probs <- matrix(NA, nrow = length(partVec), ncol = length(tauVec))
for(i in 1:length(tauVec)) {
	draws  <-rlsp(nsim, FirstPart, tauVec[i])

	parts <- apply(draws, 1, function(x) paste0(x, sep = "", collapse = ""))

	probs[, i] <- table(factor(parts, levels = partVec))/nsim 
}


res <- makeBasePartitionPlot(partitionMatrix = partMat, 
			probabilities = probs,
			c0 = FirstPart,
			name = "test LSP", 
			xLabel = "$\\tau$", xVec = tauVec)

plotRes <-  res$plot +  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Location-Scale Partition Distribution")


plotRes <- plotRes + scale_y_continuous(breaks= res$yTicks, limits = c(0, 1), expand = c(0, 0),
	                   sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
	                                       labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) 
	
pt <- plotRes + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[1] + 0.02,ymax = res$yTicks[1] + 0.02,  xmin = -0.5,xmax = 0.5) +
annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[2] + 0.02, ymax = res$yTicks[2] + 0.02,  xmin = -0.5,xmax = 0.5) +
annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[3] + 0.02,ymax = res$yTicks[3] + 0.02, xmin = -0.5,xmax = 0.5) +
annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[4] + 0.02 ,ymax = res$yTicks[4] + 0.02,  xmin = -0.5,xmax = 0.5) +
annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[5] + 0.02,ymax = res$yTicks[5] + 0.02,  xmin = -0.5,xmax = 0.5)


gt <- ggplot_gtable(ggplot_build(pt))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid_p = grid.arrange(gt)


# plotRes
ggsave(grid_p, file = paste0("figures/fig1_LSP_.pdf"), 
		device = "pdf", height = 20, width = 22, unit = "cm", dpi = 300) 

###################
## Heatmap

partMat <- dist_from(FirstPart, return_partitions = TRUE)$partitions
partMat <- partMat + 1
partVec <- apply(partMat, 1, function(x)  paste0(x, collapse = ""))

colnames(probs) <- tauVec
rownames(probs) <- partVec

probs2 <- probs[-which(rownames(probs) == "11222"), ]

probsDF <- melt(probs2)
probsDF$Var1 <- factor(probsDF$Var1)
probsDF$Var2 <- factor(probsDF$Var2)

# Color Brewer palette
ggplot(probsDF, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
    xlab("partition") + ylab("") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file = paste0("figures/fig2_LSP_heatmap.pdf"), 
    device = "pdf", height = 8, width = 32, unit = "cm", dpi = 300) 

