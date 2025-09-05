### Plots results from simulations from temporal partition process ###
### Pairwise distance between partitions ####
# library(MixSim)
library(fields)
library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
library(mclust)
theme_set(theme_minimal())
############################

## Wrapping passages in a function
plot_results <- function(res, title = "") { 

outDir <- paste0("figures/", res$name, "/")
dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)

nsim <- length(res$ARIout)

ARIoutMat <- Reduce('+', res$ARIout)/nsim
VIoutMat <- Reduce('+', res$VIout)/nsim
F1outMat <- Reduce('+', res$F1out)/nsim

 
dfARI <- reshape2::melt(ARIoutMat)
dfVI  <- reshape2::melt(VIoutMat) 
dfF1  <- reshape2::melt(F1outMat) 

# ############################
### distances over time ###
############################

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
  scale_fill_gradientn(colors= tim.colors(),limits = c(0, 1)) + 
     xlab('')+  ylab('') + 
     theme_classic(14)+  
     ggtitle("Pairwise allocation probabilities\nmarginal over time") + 
     theme( aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, linewidth=0.5), 
            legend.title = element_blank(), 
            legend.position="right")
pairPlotM
ggsave(pairPlotM, filename = paste0(outDir, "pairWiseMarginal.pdf"),
        width = 6, height = 6)

## overtime
dfOverTime  <- reshape2::melt(pairOverTime) 
dfOverTime$Var3 = dfOverTime$Var3 -1 

pairPlotOT <- ggplot(dfOverTime) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
     facet_wrap(.~as.factor(Var3), nrow = 2) + 
  scale_fill_gradientn(colors= tim.colors(),limits = c(0, 1)) + 
#  scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
#        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
      xlab('')+  ylab('') + 
     theme_classic(14)+  
     ggtitle(title) + 
     theme(axis.text.x = element_blank() , 
            axis.text.y = element_blank() , 
            aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, linewidth=0.5), 
            legend.title = element_blank(), 
            legend.position="right", 
           axis.ticks = element_blank(), 
           axis.line  =element_blank()) +
  theme(legend.key.height = unit(2, "cm"))

# pairPlotOT
ggsave(pairPlotOT, filename = paste0(outDir, "pairWiseOverTime.pdf"),
        width = 11, height = 5)


pairPlotLong <- ggplot(dfOverTime) + 
     geom_tile(aes(x = as.factor(Var1), y = as.factor(Var2), 
        fill = value), color = "white") +
     facet_wrap(.~as.factor(Var3), nrow = 1) + 
  scale_fill_gradientn(colors= tim.colors(),limits = c(0, 1)) + 
#  scale_fill_gradient2(low = "#e5f5f9",  mid = "#99d8c9", 
#        high = "#2ca25f", midpoint = 0.5, limits = c(0, 1)) + 
      xlab('')+  ylab('') + 
     theme_classic(14)+  
     theme(axis.text.x = element_blank() , 
            axis.text.y = element_blank() , 
            aspect.ratio=0.9/1, 
            panel.border = element_rect(colour = "white", fill=NA, linewidth=0.1), 
            legend.title = element_blank(), 
            legend.position="right", 
           axis.ticks = element_blank(), 
           axis.line  =element_blank()) +
  theme(legend.key.height = unit(2, "cm"))

# pairPlotOT
ggsave(pairPlotLong, filename = paste0(outDir, "pairWiseOverTimeLong.pdf"),
        width = 11, height = 5)

saveRDS(pairPlotLong, file = paste0(outDir, "pairWiseOverTime.rds"))
} 

############################
## read results
############################
res <- readRDS("results/samples_diffAlpha_constantOverTime.rds")

plot_results(res, title = "Markovian model - unit local prior")

res <- readRDS("results/samples_differentEta_phi05_AR.rds")

plot_results(res, title = "Latent autoregressive - stationary model")


res <- readRDS("results/samples_differentEta_phi08_AR.rds")

plot_results(res, title = "Latent autoregressive - stationary model")

res <- readRDS("results/samples_differentEta_phi15_AR.rds")

plot_results(res, title = "Latent autoregressive - nonstationary, positive")

res <- readRDS("results/samples_differentEta_phineg15_AR.rds")

plot_results(res, title = "Latent autoregressive - nonstationary, negative")

###
library(patchwork)
library(cowplot)

p1 <- readRDS("figures/diffAlpha_constantOverTime/pairWiseOverTime.rds")
p2 <- readRDS("figures/differentEta_phi08_AR/pairWiseOverTime.rds")
p3 <- readRDS("figures/differentEta_phi15_AR/pairWiseOverTime.rds")

prow <- plot_grid(
  p1 + theme(legend.position="none") + theme(plot.margin = unit(c(0, 0, -1, 0), "cm")) + 
  ggtitle("a) Markovian model - unit local prior"),
  p2 + theme(legend.position="none") + theme(plot.margin = unit(c(-1, 0, 0, 0), "cm")) +
  ggtitle("b) Latent autoregressive - stationary"),
  p3 + theme(legend.position="none") + theme(plot.margin = unit(c(-2, 0, 0, 0), "cm")) +
  ggtitle("c) Latent autoregressive - nonstationary"),
#  align = 'vh',
#  labels = c("(a)", "(b)", "(c)"),
  nrow = 3
)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
prow <- plot_grid(prow, legend, rel_widths = c(3, .4))
prow
# all plots
ggsave(prow, filename = paste0("figures/", "fig1_pairWiseOverTime.pdf"),
        width = 13, height = 7)
