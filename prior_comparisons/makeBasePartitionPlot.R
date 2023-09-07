makeBasePartitionPlot <- function(partitionMatrix, probabilities, c0, name, highlightColor = "#00B0DA", xLabel = "", minimal = FALSE, xVec = alphaVec){
	require(ggplot2)
	require(latex2exp)

	nBlocks = apply(partitionMatrix, 1, function(x) length(unique(x)))
	biggestBlock = apply(partitionMatrix, 1, function(x) max(table(x)))
	indBiggestBlock = order(as.numeric(paste(nBlocks, biggestBlock, sep = "")), decreasing = TRUE)
	# 
	biggestBlock[indBiggestBlock]

	tmp = probabilities[indBiggestBlock, ]
	df = reshape2::melt(tmp) ## data frame Var1 = partition index Var2 psi index 

	## Find position of the c0 partition
	c0Pos = which(apply(partitionMatrix[indBiggestBlock, ], 1, function(x) all(x == c0)))

	df$isC0 = as.numeric(df$Var1==c0Pos)

	g <- grouping(nBlocks[indBiggestBlock])

	yTicks = 1 - c(cumsum(tmp[,1])[g[attr(g, "ends")]], 0)
	p = ggplot(df, aes(x = Var2, y = value)) + 
	geom_area(aes(group = Var1, fill = factor(isC0)),col = "grey20" ) +
	scale_fill_manual(values = c('grey90',highlightColor))+
	theme_bw(18) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
	scale_x_continuous(breaks = 1:length(xVec),labels = xVec, expand = c(0, 0))+
		theme(legend.position = 'none',  plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
		#theme(plot.margin=grid::unit(c(1,1,1,1), "cm"))
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		  axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
	      axis.ticks.y.left = element_line(linewidth = 0.6, colour = "grey20", linetype = "dotted"),
	      axis.ticks.length.y.left  = unit(2, "cm"),
	      axis.ticks.length.x.bottom  = unit(2, "cm"),
	      axis.text.y.right = element_text( size = rel(1.6)),
	      axis.ticks.x = element_blank(),
	      axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
	      axis.title.x=element_text(size=20,face="bold"),
	       axis.title.y.right=element_text(size=20)
	      ) + xlab(TeX(xLabel))	

	out <-list(yTicks = yTicks, plot = p)
}
