################################
## this script plots results from simulation
################################
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(latex2exp)

plotColors <- brewer.pal(3, "Dark2")

resList <- list()
for(i in 1:9) resList[[i]] <- readRDS(paste0("output/scenario_", i, ".rds"))

dfIP1 <- data.frame()

for(i in 1:9) {
	dfIP <- as.data.frame(resList[[i]]$ariTruePartIP)
	dfIP$meanSep <- paste0("mean separation = ", resList[[i]]$scenario$meanSep)
	dfIP$initialPartition <- paste0("initial partition = ", resList[[i]]$scenario$initialPartition)
	dfIPlong <- melt(dfIP, id.vars = c("initialPartition", "meanSep"))
	dfIP1 <-rbind(dfIP1, dfIPlong)
}

# labels <- c(TeX("initial partition = $\\rho_{merge}$" ), 
# TeX("initial partition = $\\rho_{split}$" ), 
# TeX("initial partition = $\\rho_{truth}$" ))

# dfIP1$initialPartition <- factor(dfIP1$initialPartition, labels = labels)

p1 <- ggplot(dfIP1, aes(y = value, x = variable, color = meanSep)) + geom_boxplot() + theme_bw() +
	facet_grid(.~ initialPartition) + 
#	scale_colour_manual(values = plotColors) +
	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\alpha")) + ggtitle("informed partition") + theme(legend.position = "none")
p1



ggsave(filename = "figXX_simulationIP.pdf", plot = p1, width = 11, height = 4)

dfCP1 <- data.frame()

for(i in 1:9) {
	dfCP <- as.data.frame(resList[[i]]$ariTruePartCP)
	dfCP$meanSep <- paste0("mean separation = ", resList[[i]]$scenario$meanSep)
	dfCP$initialPartition <- paste0("initial partition = ", resList[[i]]$scenario$initialPartition)
	dfCPlong <- melt(dfCP, id.vars = c("initialPartition", "meanSep"))
	dfCP1 <-rbind(dfCP1, dfCPlong)
}

p2<- ggplot(dfCP1, aes(y = value, x = variable,  color = meanSep)) + geom_boxplot() + theme_bw() +
	facet_grid(.~ initialPartition ) + 
	ylab(TeX("ARI($\\hat{\\rho}, \\rho_{0}$)")) + xlab(TeX("\\psi")) + ggtitle("centered partition process") + theme(legend.position = "bottom")

ggsave(filename = "figXX_simulationCP.pdf", plot = p2, width = 11, height = 4.5) 


# resList[[1]]$ariInitialPartCP

# plot(res$ariInitialPartIP[, 1], res$ariInitialPartCP[, 1])