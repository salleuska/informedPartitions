####################
## Sally Paganin
## 2023/05/19
####################
## This script produces a density plot over the set partition space of 5 elements
## for the centered partition process 
####################
library(CPLogit)
library(ggplot2)
library(grid)
library(reshape2)
library(gridExtra)
####################
# Bell numbers - sequence A000110
bell_n <- c(1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597, 27644437, 190899322, 1382958545, 10480142147, 82864869804)
####################
## Functions for the CP prior
####################
CPpenalization <- function(distance, psi, returnLog = FALSE) {
  if(returnLog){
    -psi*distance
  } else {
    exp(-psi*distance)
  }
}

####################
## Dirichlet Process EPPF e
####################
dpEPPF <- function(partition, concentration, returnLog = FALSE){
  nObs  <- length(partition)
  nBlocks  <- length(unique(partition))
  blockSizes  <- table(partition)
  
  logProb  <- nBlocks*log(concentration) + sum(lgamma(blockSizes))
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}

###########################
## Finite Dirichlet EPPF
###########################

finiteDirichletEPPF <- function(partition, concentration, nClusters, returnLog = FALSE){
  nObs  <- length(partition)
  nBlocks  <- length(unique(partition))
  blockSizes  <- table(partition)
  
  log_c  <- 0
  if(nBlocks < nClusters) {log_c  <- - lgamma(nClusters - nBlocks +1)}
  logProb  <-   lgamma(nClusters +1) + log_c + sum(lgamma(concentration/nClusters + blockSizes) - lgamma(concentration/H))
  ##  + lgamma(concentration) - lgamma(concentration + nObs) # (normalizaoitn constant)
  
  if(returnLog) {
     logProb   
  } else {
       exp(logProb)  
  }
}

##############################
## Pitman-Yor process eppf
##############################
pyEPPF <- function(partition, concentration, discount, returnLog = FALSE){
  n  <-length(partition)
  nBlocks  <-length(unique(partition))
  blockSizes  <-table(partition)
  
  if(nBlocks == 1){
    logProb  <--
      lgamma(n + concentration) + lgamma(concentration + 1) + 
      sum(lgamma(discount + blockSizes) - lgamma(1 - discount))   

  } else {
    logProb  <- sum(log(concentration + 1:(nBlocks-1)*discount)) - 
      lgamma(n + concentration) + lgamma(concentration + 1) + 
      sum(lgamma(discount + blockSizes) - lgamma(1 - discount))   
  }
  
  
  if(returnLog) {
     logProb   
  } else {
     exp(logProb)  
  }
}

##################
## Computation and plotting
##################
## Settings ----
## Set partition color for c0
col2 <- "#00B0DA"
alpha <- 1  ## DP concentration parameter
# psiVal <- c(0,0.5,1,1.5,2,2.5,3, 3.5, 4, 4.5, 5)
psiVal <- seq(0, 10,length = 11)
#psiVal <- seq(0, 15,length = 16)

c0 <- c(1,1,2,2,2) -1 
####################
## Computation -----
## Generate all possible partitions and compute distances from c0 and g
dd <- CPLogit::dist_from(c0, return_partitions = TRUE)
## Numerator of CP process penalization
CPNum <- sapply(psiVal, function(x) CPpenalization(dd$distances, psi = x))
## Numerator of DP process EPP
DPNum  <-apply(dd$partitions, 1 ,function(x) dpEPPF(partition = x, concentration =  alpha))
## Compute normalized partition probabilities
partitionProbs  <-apply(CPNum, 2, function(x) (DPNum*x)/sum(DPNum*x))

####################
## Plotting ----

### sort partition by block size
nBlocks  <-apply(dd$partition, 1, function(x) length(unique(x)))
biggestBlock  <-apply(dd$partition, 1, function(x) max(table(x)))
indBiggestBlock <- order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)


tmp  <-partitionProbs[indBiggestBlock, ]
df  <-melt(tmp) ## data frame Var1 = partition index Var2 psi index 

## Find position of the c0 partition
c0Pos  <-which(apply(dd$partitions[indBiggestBlock, ], 1, function(x) all(x == c0)))
df$isC0  <-as.numeric(df$Var1==c0Pos)

df$name <- "Centered Partition Process"

g <- grouping(nBlocks[indBiggestBlock])

yTicks <- 1 - c(cumsum(tmp[,1])[g[attr(g, "ends")]], 0)
# p <- ggplot(df, aes(x = Var2, y = value)) + 
#   geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
#   scale_fill_manual(values = c('grey90',col2))+
#   theme_minimal(18) +
#   scale_x_continuous(breaks = 1:length(psiVal),labels = psiVal) +
#   theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
#   scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
#                      sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
#                                          labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
#   theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
#         axis.ticks.y.left = element_line(linewidth = 0.6, colour = "grey20", linetype = "dotted"),
#         axis.ticks.length.y.left  = unit(2, "cm"),
#         axis.ticks.length.x.bottom  = unit(2, "cm"),
#         axis.text.y.right = element_text( size = rel(1.6)),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
#         axis.title.x=element_text(size=20,face="bold"),
#          axis.title.y.right=element_text(size=20),
#          plot.title = element_text(hjust = 0.5)
#         ) +
#   xlab(expression(psi)) + ggtitle("Centered Partition Process")
# p

p <- ggplot(df, aes(x = Var2, y = value)) + 
  geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
  scale_fill_manual(values = c('grey90',col2))+
  theme_bw(18) +
  scale_x_continuous(breaks = 1:length(psiVal),labels = psiVal, expand = c(0, 0)) +
  theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
  scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
                     sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                         labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
        axis.ticks.y.left = element_line(linewidth = 0.6, colour = "grey20", linetype = "dotted"),
        axis.ticks.length.y.left  = unit(2, "cm"),
        axis.ticks.length.x.bottom  = unit(2, "cm"),
        axis.text.y.right = element_text( size = rel(1.6)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
        axis.title.x=element_text(size=20,face="bold"),
         axis.title.y.right=element_text(size=20),
         plot.title = element_text(hjust = 0.5)
        ) +
  xlab(expression(psi)) + ggtitle("Centered Partition Process")
#p

pt <- p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = yTicks[1] + 0.02,ymax = yTicks[1] + 0.02,  xmin = -0.6,xmax = 0.6) +
annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = yTicks[2] + 0.02, ymax = yTicks[2] + 0.02,  xmin = -0.6,xmax = 0.6) +
annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = yTicks[3] + 0.02,ymax = yTicks[3] + 0.02, xmin = -0.6,xmax = 0.6) +
annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = yTicks[4] + 0.02 ,ymax = yTicks[4] + 0.02,  xmin = -0.6,xmax = 0.6) +
annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = yTicks[5] + 0.02,ymax = yTicks[5] + 0.02,  xmin = -0.6,xmax = 0.6)

gt <- ggplot_gtable(ggplot_build(pt))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid_p = grid.arrange(gt)

ggsave(grid_p, file = paste0("figures/fig1_CP_priorDensity.pdf"), 
    device = "pdf", height = 20, width = 22, unit = "cm", dpi = 300) 

####################################################################################################
## Heatmap other partitions

partMat <- dist_from(c0, return_partitions = TRUE)$partitions
partMat <- partMat + 1
partVec <- apply(partMat, 1, function(x)  paste0(x, collapse = ""))

colnames(partitionProbs) <- psiVal
rownames(partitionProbs) <- partVec

saveRDS(partitionProbs, "CPP_partitionProbs.rds")
probs <- partitionProbs[-which(rownames(partitionProbs) == "11222"), ]

probsDF <- melt(probs)
probsDF$Var1 <- factor(probsDF$Var1)
probsDF$Var2 <- factor(probsDF$Var2)

# Color Brewer palette
ggplot(probsDF, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
    xlab("partition") + ylab("") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file = paste0("figures/fig2_CP_heatmap.pdf"), 
    device = "pdf", height = 8, width = 32, unit = "cm", dpi = 300) 



