################################
## this script plots results from posterior simulation
################################
library(reshape2)
library(ggplot2)
library(latex2exp)
library(dplyr)
library(purrr)
## library(mclust)   ## adjustedRandIndex - not loaded because it masks purrr::map
library(patchwork)
################################
plotColors <- c("#F8766D", "#54B502", "#08B9EC")
################################
## create a list with all the results
resList <- list()
for(i in 1:9) resList[[i]] <- readRDS(paste0("output/scenario_", i, ".rds"))
################################
## function to extract results and collect them in data frame with long format for plotting

extractRes <- function(dataList, measureName) { 
	tmp <-  purrr::map(resList, \(x) as.data.frame(x[[measureName]]) |> 
	 							mutate(meanSep = paste0("h = ", x$scenario$meanSep), 
	 								   initialPartition = x$scenario$initialPartition)) |> bind_rows()
	
	tmplong <- melt(tmp, id.vars = c("initialPartition", "meanSep"))

	tmplong$initialPartition<- as.factor(tmplong$initialPartition)
	levels(tmplong$initialPartition) <- c(TeX("$\\rho_{0} = \\rho_{merge}$"), TeX("$\\rho_{0} = \\rho_{split}$"), TeX("$\\rho_{0} = \\rho_{true}$"))
	tmplong$initialPartition<- factor(tmplong$initialPartition, levels = levels(tmplong$initialPartition)[c(3,1,2)] )
	tmplong
}

################################
## Create a dataframe containing results for each prior - ari from true partition
dfIP  <- extractRes(resList, "ariTruePartIP")
dfCP  <- extractRes(resList, "ariTruePartCP")
dfLSP <- extractRes(resList, "ariTruePartLSP")
################################
## info about initial partitions
N <- 100
trueP  <-  rep(1:4, each=N/4) 
mergeP <-  rep(1:2, each=N/2)
splitP <-  c(rep(1, 12), rep(2, 13), rep(3, 12), rep(4, 13),
	rep(5, 12), rep(6, 13), rep(7, 12), rep(8, 13))

distFromTruePart <- data.frame(initialPartition = c("true", "merge", "split"), 
							   dist = c(mclust::adjustedRandIndex(trueP, trueP), mclust::adjustedRandIndex(trueP, mergeP), mclust::adjustedRandIndex(trueP, splitP)))

distFromTruePart$initialPartition <- factor(distFromTruePart$initialPartition, labels = c(TeX("$\\rho_{0} = \\rho_{merge}$"), TeX("$\\rho_{0} = \\rho_{split}$"), TeX("$\\rho_{0} = \\rho_{true}$")))
################################


plotIP <- ggplot(dfIP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	geom_hline(data =distFromTruePart, aes(yintercept = dist), linetype = "dashed") + 
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	ylab(TeX("E(ARI($\\rho, \\rho_{true}$)| Y)")) + xlab(TeX("\\alpha")) + 
	ggtitle("Informed Partition Model") + theme(legend.position = "none") 

plotIP 

plotCP <- ggplot(dfCP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	geom_hline(data =distFromTruePart, aes(yintercept = dist), linetype = "dashed") + 
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	ylab(TeX("E(ARI($\\rho, \\rho_{true}$)| Y)"))  + xlab(TeX("\\psi")) + 
	ggtitle("Centered Partition Process") + theme(legend.position = "none")

plotCP 


labelsLSP <- c(c(10, 5, 1, 0.05), TeX("\\frac{1}{mlog(m)}"), TeX("\\frac{0.1}{mlog(m)}"))

plotLSP <- ggplot(dfLSP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	geom_hline(data =distFromTruePart, aes(yintercept = dist), linetype = "dashed") + 
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	scale_x_discrete(labels= labelsLSP) + 
	ylab(TeX("E(ARI($\\rho, \\rho_{true}$)| Y)"))  + xlab(TeX("\\nu")) + 
	ggtitle("Location-Scale Partition Distribution") +
	theme(legend.position = "bottom", legend.title = element_blank(),  
		axis.text.x = element_text(angle = 0))


p <- plotIP / plotCP / plotLSP 

ggsave(filename = "fig_ARI_PostComparison.pdf", plot = p, width = 11, height = 12.5) 



# #############
## Create a dataframe containing results for each prior - 
dfIPWaic  <- extractRes(resList, "waicIP")
dfCPWaic  <- extractRes(resList, "waicCP")

dfIPlpml  <- extractRes(resList, "lmplIP")
dfCPlpml  <- extractRes(resList, "lmplCP")
dfLSPlpml <- extractRes(resList, "lmplLSP")

################################
plotIPlpml <- ggplot(dfIPlpml, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	ylab("LPML") + xlab(TeX("\\alpha")) + 
	ggtitle("Informed Partition Model") + theme(legend.position = "none") 

plotCPlpml<- ggplot(dfCPlpml, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	ylab("LPML") + xlab(TeX("\\alpha")) + 
	ggtitle("Centered Partition Process") + theme(legend.position = "none") 

plotLSPlpml <- ggplot(dfLSPlpml, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	scale_x_discrete(labels= labelsLSP) + 
	ylab(TeX("LPML"))  + xlab(TeX("\\nu")) + 
	ggtitle("Location-Scale Partition Distribution") +
	theme(legend.position = "bottom", legend.title = element_blank(),  
		axis.text.x = element_text(angle = 0))

# plpml <- plotIPlpml / plotCPlpml / plotLSPlpml 
# ggsave(filename = "fig_LPML_postComp.pdf", plot = plpml, width = 11, height = 12.5) 

# ###

# plotIPWaic <- ggplot(dfIPWaic, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	scale_colour_manual(values = plotColors)  +
# 	ylab("WAIC") + xlab(TeX("\\alpha")) + 
# 	ggtitle("Informed Partition Model") + theme(legend.position = "none") 


# plotCPWaic<- ggplot(dfCPWaic, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	scale_colour_manual(values = plotColors)  +
# 	ylab("WAIC") + xlab(TeX("\\alpha")) + 
# 	ggtitle("Centered Partition Process") + theme(legend.position = "none") 

# plotLSPWaic <- ggplot(dfLSPWaic, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition) + 
# 	scale_colour_manual(values = plotColors)  +
# 	scale_x_discrete(labels= labelsLSP) + 
# 	ylab(TeX("WAIC"))  + xlab(TeX("\\nu")) + 
# 	ggtitle("Location-Scale Partition Distribution") +
# 	theme(legend.position = "bottom", legend.title = element_blank(),  
# 		axis.text.x = element_text(angle = 0))

# pWaic <- plotIPWaic / plotCPWaic / plotLSPWaic 


# ggsave(filename = "fig_WAIC_postComp.pdf", plot = pWaic, width = 11, height = 12.5) 
