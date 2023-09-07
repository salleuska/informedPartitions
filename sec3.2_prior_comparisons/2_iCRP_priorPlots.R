#################
## Sally Paganin
## 2023/14/03 
#################
## This script produces a density plot over the set partition space of 5 elements
## for the informed partition model 
#################
library(CPLogit)
library(reshape2)
require(grid)
require(gridExtra)
source("utils/iCRP_densityFunction.R")
source("utils/makeBasePartitionPlot.R")
####################################
## create figures and output directory if missing
dir.create(file.path("output"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures"), recursive = TRUE, showWarnings = FALSE)
####################################

## values of alpha parameter
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

res <- makeBasePartitionPlot(partitionMatrix = partMat, 
			probabilities = t(probabilities),
			c0 = FirstPart,
			name = name,
			xLabel = "$\\alpha$")

plotRes <-  res$plot +  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("iCRP")


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
ggsave(grid_p, file = paste0("figures/fig1_iCRP_", name, ".pdf"), 
		device = "pdf", height = 20, width = 22, unit = "cm", dpi = 300) 

####################################
### HEATMAP

partMat <- dist_from(FirstPart, return_partitions = TRUE)$partitions
partMat <- partMat + 1
partVec <- apply(partMat, 1, function(x)  paste0(x, collapse = ""))


colnames(probabilities) <- partVec
rownames(probabilities) <- alphaMat[,1]

saveRDS(probabilities, file = "output/iCRP_partitionProbs.rds")

probs <- t(probabilities[-which(colnames(probabilities) == "11222"), ])

probsDF <- melt(probs)
probsDF$Var1 <- factor(probsDF$Var1)
probsDF$Var2 <- factor(probsDF$Var2)

# Color Brewer palette
ggplot(probsDF, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
    xlab("partition") + ylab("") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file = paste0("figures/fig2_iCRP_heatmap.pdf"), 
    device = "pdf", height = 8, width = 32, unit = "cm", dpi = 300) 



####################################
## Example 2 - changing probabilities for bigger cluster
####################################
################################
## create an empty list
fixedAlphaVals <- c(0, 0.25, 0.5, 0.75, 1)

## settings
# cluster 1 probs = 0
FirstPart<- c(1,1,2,2,2) 
## Define scenario for the pegging probabilities
nVals <- 11
alphaMat <- matrix(NA, nrow = nVals, ncol = length(FirstPart))
alphaMat[, 3:5] <- seq(0, 1, length = nVals)

computedProbs <- list()
for(i in 1:length(fixedAlphaVals)) {
	alphaMat[, 1:2] <- fixedAlphaVals[i]

	computedProbs[[i]] <- apply(partMat, 1,  function(x) apply(alphaMat, 1, function(y) dIDP(x, part0 = FirstPart, concentration = M, alpha = y)))

}
names(computedProbs) <- fixedAlphaVals

## saving to avoid re-computing
# saveRDS(computedProbs, file = "output/iCRP_unitAlphaProbs.rds")
# computedProbs	<- readRDS("output/iCRP_unitAlphaProbs.rds")
########
## Settings 
c0 = FirstPart
highlightColor = "#00B0DA"
xLabel = "$\\alpha$"
## common to all plots
nBlocks = apply(partMat, 1, function(x) length(unique(x)))
biggestBlock = apply(partMat, 1, function(x) max(table(x)))
indBiggestBlock = order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)
biggestBlock[indBiggestBlock]
##########################
plotList <- list()
for(i in 1:length(fixedAlphaVals)){ 
	require(ggplot2)
	require(latex2exp)
	probabilities <- t(computedProbs[[i]])
	tmp = probabilities[indBiggestBlock, ]
	df = reshape2::melt(tmp) ## data frame Var1 = partition index Var2 psi index 

	## Find position of the c0 partition
	c0Pos = which(apply(partMat[indBiggestBlock, ], 1, function(x) all(x == c0)))

	df$isC0 = as.numeric(df$Var1==c0Pos)

	g <- grouping(nBlocks[indBiggestBlock])

	
	p <- ggplot(df, aes(x = Var2, y = value)) + 
	geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
	scale_fill_manual(values = c('grey90',highlightColor))+
	theme_bw(18) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
	scale_x_continuous(breaks = 1:length(seq(0, 1, length = nVals)),labels = seq(0, 1, length = nVals), expand = c(0, 0))+
	  #	scale_x_continuous(breaks = 1:11, labels = c(0, rep("", 4), 0.5, rep("", 4), 1),  expand = c(0, 0))+
		theme(legend.position = 'none',  plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
		#theme(plot.margin=grid::unit(c(1,1,1,1), "cm"))
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		  axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
		  axis.ticks.y.left = element_line(linewidth = 0.6, colour = "grey20", linetype = "dotted"),
		  axis.ticks.length.y.left  = unit(2, "cm"),
	      axis.ticks.length.x.bottom  = unit(2, "cm"),
	      axis.text.y.right = element_text( size = rel(1.6)),
	      axis.ticks.x = element_blank(),
	      axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1)),
	      axis.title.x=element_text(size=20,face="bold"),
	       axis.title.y.right=element_text(size=20)
	      ) + xlab(TeX(xLabel))	
	
	  title <- paste0("iCRP with $\\alpha_{1} = \\alpha_{2} = ", fixedAlphaVals[i], "\\; \\alpha_{3} = \\alpha_{4} = \\alpha_{5} = \\alpha$")
	
  p  <- p+  theme(plot.title = element_text(hjust = 0.5)) + ggtitle(TeX(title))
	
	  p <- p + scale_y_continuous(breaks= res$yTicks, limits = c(0, 1), expand = c(0, 0),
	                                        sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
	                                                            labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) 
	
	  p <- p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[1] + 0.02,ymax = res$yTicks[1] + 0.02,  xmin = -0.5,xmax = 0.5) +
	  annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[2] + 0.02, ymax = res$yTicks[2] + 0.02,  xmin = -0.5,xmax = 0.5) +
	  annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[3] + 0.02,ymax = res$yTicks[3] + 0.02, xmin = -0.5,xmax = 0.5) +
	  annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[4] + 0.02 ,ymax = res$yTicks[4] + 0.02,  xmin = -0.5,xmax = 0.5) +
	  annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = res$yTicks[5] + 0.02,ymax = res$yTicks[5] + 0.02,  xmin = -0.5,xmax = 0.5)
	
	
	gt <- ggplot_gtable(ggplot_build(p))
	gt$layout$clip[gt$layout$name == "panel"] <- "off"
	grid_p = grid.arrange(gt)

	
	# plotRes
	ggsave(grid_p, file = paste0("figures/figApp_iCRP_fixed_", fixedAlphaVals[i], ".pdf"), 
	       device = "pdf", height = 20, width = 22, unit = "cm", dpi = 300) 
	
}
# 	
# library(cowplot)
# plotToSave <-  plot_grid(plotlist = plotList, nrow = 1)
# 
# ggsave(plotToSave, file = paste0("figures/iCRP_partial.pdf"), 
# 		device = "pdf", height = 15, width = 75, unit = "cm", dpi = 300) 
# 
