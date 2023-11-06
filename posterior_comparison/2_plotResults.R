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
	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{true}$)")) + xlab(TeX("\\alpha")) + 
	ggtitle("Informed Partition Model") + theme(legend.position = "none") 

plotIP 

plotCP <- ggplot(dfCP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	geom_hline(data =distFromTruePart, aes(yintercept = dist), linetype = "dashed") + 
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{true}$)")) + xlab(TeX("\\psi")) + 
	ggtitle("Centered Partition Process") + theme(legend.position = "none")

plotCP 


labelsLSP <- c(c(10, 5, 1, 0.05), TeX("\\frac{1}{mlog(m)}"), TeX("\\frac{0.1}{mlog(m)}"))

plotLSP <- ggplot(dfLSP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	geom_hline(data =distFromTruePart, aes(yintercept = dist), linetype = "dashed") + 
	facet_grid(.~ initialPartition, labeller=label_parsed) + 
	scale_colour_manual(values = plotColors)  +
	scale_x_discrete(labels= labelsLSP) + 
	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{true}$)")) + xlab(TeX("\\tau")) + 
	ggtitle("Location-Scale Partition Distribution") +
	theme(legend.position = "bottom", legend.title = element_blank(),  
		axis.text.x = element_text(angle = 0))


p <- plotIP / plotCP / plotLSP 

ggsave(filename = "figXX_ARI_PostComparison.pdf", plot = p, width = 11, height = 12.5) 



# ## Create a dataframe containing results for each prior - ari from true partition
# dfIP  <- extractRes(resList, "ariInitialPartIP")
# dfCP  <- extractRes(resList, "ariInitialPartCP")
# dfLSP <- extractRes(resList, "ariInitialPartLSP")

# plotIP <- ggplot(dfIP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	scale_colour_manual(values = plotColors)  +
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\alpha")) + 
# 	ggtitle("Informed Partition Model") + theme(legend.position = "none") 


# plotCP <- ggplot(dfCP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	scale_colour_manual(values = plotColors)  +
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\psi")) + 
# 	ggtitle("Centered Partition Process") + theme(legend.position = "none")


# labelsLSP <- c(c(10, 5, 1, 0.05), TeX("\\frac{1}{mlog(m)}"), TeX("\\frac{0.1}{mlog(m)}"))

# plotLSP <- ggplot(dfLSP, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	scale_colour_manual(values = plotColors)  +
# 	scale_x_discrete(labels= labelsLSP) + 
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\tau")) + 
# 	ggtitle("Location-Scale Partition Distribution") +
# 	theme(legend.position = "bottom", legend.title = element_blank(),  
# 		axis.text.x = element_text(angle = 0))


# p <- plotIP / plotCP / plotLSP 
# ggsave(filename = "figXX_initialPArt.pdf", plot = p, width = 11, height = 12.5) 

# ###########

# for(i in 1:9) {
# 	dfIPtmp <- as.data.frame(resList[[i]]$ariTruePartIP)
# 	dfIPtmp$meanSep <- paste0("mean separation = ", resList[[i]]$scenario$meanSep)
# 	dfIPtmp$initialPartition <- paste0("initial partition = ", resList[[i]]$scenario$initialPartition)
# 	dfIPlong <- melt(dfIPtmp, id.vars = c("initialPartition", "meanSep"))
# 	dfIP1 <-rbind(dfIP1, dfIPlong)

# }

# # labels <- c(TeX("initial partition = $\\rho_{merge}$" ), 
# # TeX("initial partition = $\\rho_{split}$" ), 
# # TeX("initial partition = $\\rho_{truth}$" ))
# # dfIP1$initialPartition <- factor(dfIP1$initialPartition, labels = labels)

# p1 <- ggplot(dfIP1, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition) + 
# 	scale_colour_manual(values = c(2,3,4)) +
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\alpha")) + ggtitle("informed partition") + theme(legend.position = "none")
# p1



# ggsave(filename = "figXX_simulationIP.pdf", plot = p1, width = 11, height = 4)

# dfCP1 <- data.frame()

# for(i in 1:9) {
# 	dfCP <- as.data.frame(resList[[i]]$ariTruePartCP)
# 	dfCP$meanSep <- paste0("mean separation = ", resList[[i]]$scenario$meanSep)
# 	dfCP$initialPartition <- paste0("initial partition = ", resList[[i]]$scenario$initialPartition)
# 	dfCPlong <- melt(dfCP, id.vars = c("initialPartition", "meanSep"))
# 	dfCP1 <-rbind(dfCP1, dfCPlong)
# }

# p2<- ggplot(dfCP1, aes(y = value, x = variable,  color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition ) + 
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\psi")) + ggtitle("centered partition process") + theme(legend.position = "bottom")

# p2

# ggsave(filename = "figXX_simulationCP.pdf", plot = p2, width = 11, height = 4.5) 



# dfLSP <- data.frame()
# for(i in 1:9) {
# 	dfLSPtmp <- as.data.frame(resList[[i]]$ariTruePartLSP)
# 	dfLSPtmp$meanSep <- paste0("h = ", resList[[i]]$scenario$meanSep)
# 	dfLSPtmp$initialPartition <- paste0("initial partition = ", resList[[i]]$scenario$initialPartition)
# 	dfLSPlong <- melt(dfLSPtmp, id.vars = c("initialPartition", "meanSep"))
# 	dfLSP <-rbind(dfLSP, dfLSPlong)
# }

# dfLSP$initialPartition<- as.factor(dfLSP$initialPartition)

# levels(dfLSP$initialPartition) <- c(TeX("$\\rho_{0} = \\rho_{merge}$"), TeX("$\\rho_{0} = \\rho_{split}$"), TeX("$\\rho_{0} = \\rho_{true}$"))

# dfLSP$initialPartition<- factor(dfLSP$initialPartition, levels = levels(dfLSP$initialPartition)[c(3,1,2)] )

# p3 <- ggplot(dfLSP, aes(y = value, x = variable,  color = meanSep)) + geom_boxplot() + theme_bw() +
# 	facet_grid(.~ initialPartition, labeller=label_parsed) + 
# 	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\tau")) + ggtitle("LSP prior") + theme(legend.position = "bottom", legend.title = element_blank())
# p3


# # resList[[1]]$ariInitialPartCP

# # plot(res$ariInitialPartIP[, 1], res$ariInitialPartCP[, 1])