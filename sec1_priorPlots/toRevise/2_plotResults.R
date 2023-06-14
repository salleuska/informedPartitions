### Plots results from simulations from temporal partition process ###
### Pairwise distance between partitions ####
library(MixSim)
library(fields)
library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
theme_set(theme_minimal())
############################
## read results
############################
# res <- readRDS("results/samplesDists_firstPart_sameRho.rds")
# res <- readRDS("results/samplesDists_firstPart_diffRho_constantOverTime.rds")
#res <- readRDS("results/samplesDists_firstPart_diffRho_constantOverTime_05.rds")
res <- readRDS("results/samplesDists_firstPart_diffRho_constantOverTime_high.rds")
# res <- readRDS("results/samplesDists_firstPart_clusterSpecificRho_constantOverTime.rds")
outDir <- paste0("figures/", res$name, "/")
dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)

str(res)
nsim <- length(res$ARIout)

ARIoutMat <- Reduce('+', res$ARIout)/nsim
VIoutMat <- Reduce('+', res$VIout)/nsim
F1outMat <- Reduce('+', res$F1out)/nsim

 
dfARI <- reshape2::melt(ARIoutMat)
dfVI  <- reshape2::melt(VIoutMat) 
dfF1  <- reshape2::melt(F1outMat) 
       

ARIp <- ggplot(dfARI) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
    scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
      xlab('')+  ylab('') + 
     theme_classic(14)+  
     theme( aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, size=0.5), 
            legend.title = element_blank(), 
            legend.position="right")

VIp <- ggplot(dfVI) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = 1 - value), color = "white") +
    scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
     xlab('')+ ylab('') + 
     theme_classic(14)+  
     theme( aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, size=0.5), 
            legend.title = element_blank(), 
            legend.position="right")

F1p <- ggplot(dfF1) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
    scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
     xlab('')+
     ylab('')+
     theme_classic(14)+  
     theme( aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, size=0.5), 
            legend.title = element_blank(), 
            legend.position="right")
F1p

prow <- plot_grid(
  ARIp + theme(legend.position="none") + ggtitle('ARI') ,
  VIp + theme(legend.position="none") + ggtitle('1 - VI') ,
  F1p + theme(legend.position="none") + ggtitle('F1') ,
  align = 'vh',
  hjust = -1,
  nrow = 1
)
prow
legend <- get_legend(ARIp)
plot <- plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(plot, filename = paste0(outDir, "distanceBetweenPartition.pdf"),
        width = 10, height = 5)

# ntime <- res$settings$ntime
# par(mfrow=c(1,3))

# image.plot(ARIoutMat,  axes=FALSE, main = "ARI")
# mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
# mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

# image.plot(1 - VIoutMat,  axes=FALSE, main = "1 - VI")
# mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
# mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)

# image.plot(F1outMat,  axes=FALSE, main = "F1 distance")
# mtext(text=c(paste("",1:ntime)), side=2, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
# mtext(text=c(paste("",1:ntime)), side=1, line=0.3, at=seq(0,1,length=ntime), las=1, cex=0.8)
############################
### distances over time ###
############################

distDF <- reshape2::melt(res$VIFirstout)

distDF = distDF |>
	mutate(time = factor(Var1, levels =  1:length(unique(distDF$Var1))))

distVI <- ggplot(filter(distDF, time != "1")) +
    geom_density_ridges(aes(x = value, y = time),
    stat = "binline", bins = 100, scale = 0.8,
    draw_baseline = FALSE) + 
    xlim(c(0, 1)) + 
    xlab("variation of information") + 
    ggtitle("Distance from baseline partition (T = 1)")

ggsave(distVI, filename = paste0(outDir,"distFromFirstVI.pdf"),
        width = 6, height = 8)

##############
distDF <- reshape2::melt(res$F1Firstout)

distDF = distDF |>
	mutate(time = factor(Var1, levels =  1:length(unique(distDF$Var1))))

distF1 <- ggplot(filter(distDF, time != "1")) +
    geom_density_ridges(aes(x = 1 - value, y = time),
    stat = "binline", bins = 100, scale = 0.8,
    draw_baseline = FALSE) + 
    xlim(c(0, 1)) + 
    xlab("1 - F1 distance") + 
    ggtitle("Distance from baseline partition (T = 1)")

ggsave(distF1, filename = paste0(outDir,"distFromFirstF1.pdf"),
        width = 6, height = 8)

##############
distDF <- reshape2::melt(res$ARIFirstout)

distDF = distDF |>
	mutate(time = factor(Var1, levels =  1:length(unique(distDF$Var1))))

distARI <- ggplot(filter(distDF, time != "1")) +
    geom_density_ridges(aes(x =  value, y = time),
    stat = "binline", bins = 100, scale = 0.8,
    draw_baseline = FALSE) + 
    xlab("ARI distance") + 
    ggtitle("Distance from baseline partition (T = 1)")

ggsave(distARI, filename = paste0(outDir,"distFromFirstARI.pdf"),
        width = 6, height = 8)
###################
## Pariwise allocations over time
###################

dim(res$pairWise)
pairOverTime <- apply(res$pairWise, c(1,2,3), mean)
pairMarginal <- apply(res$pairWise, c(1,2), mean)

dfMarginal  <- reshape2::melt(pairMarginal) 
       
pairPlotM <- ggplot(dfMarginal) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
    scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
      xlab('')+  ylab('') + 
     theme_classic(14)+  
     ggtitle("Pairwise allocation probabilities\nmarginal over time") + 
     theme( aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, size=0.5), 
            legend.title = element_blank(), 
            legend.position="right")
pairPlotM
ggsave(pairPlotM, filename = paste0(outDir, "pairWiseMarginal.pdf"),
        width = 6, height = 6)


dfOverTime  <- reshape2::melt(pairOverTime) 
pairPlotOT <- ggplot(dfOverTime) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
     facet_wrap(.~as.factor(Var3), nrow = 2) + 
    scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
      xlab('')+  ylab('') + 
     theme_classic(14)+  
     ggtitle("Pairwise allocation probabilities - over time") + 
     theme(axis.text.x = element_blank() , 
            axis.text.y = element_blank() , 
            aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, size=0.5), 
            legend.title = element_blank(), 
            legend.position="right")
pairPlotOT
ggsave(pairPlotOT, filename = paste0(outDir, "pairWiseOverTime.pdf"),
        width = 12, height = 7)
