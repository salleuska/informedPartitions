makePartitionPlot <- function(partitionMatrix, probabilities, c0, name, highlightColor = "#00B0DA", xLabel = ""){
	require(ggplot2)
	require(grid)
	require(gridExtra)
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
	theme_minimal(18) +
	scale_x_continuous(breaks = 1:length(alphaVec),labels = alphaVec)+
		theme(legend.position = 'none',plot.margin=grid::unit(c(1,1,1,1), "cm")) + 
	scale_y_continuous(breaks= yTicks, limits = c(0, 1), expand = c(0, 0),
	                   sec.axis = sec_axis(trans = ~., name = "Cumulative probability", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
	                                       labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"))) +
	theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), 
	      axis.ticks.y.left = element_line(size = 0.6, colour = "grey20", linetype = "dotted"),
	      axis.ticks.length.y.left  = unit(2, "cm"),
	      axis.ticks.length.x.bottom  = unit(2, "cm"),
	      axis.text.y.right = element_text( size = rel(1.6)),
	      axis.ticks.x = element_blank(),
	      axis.text.x = element_text(margin = margin(t = -1.5, unit = "cm"), size = rel(1.6)),
	      axis.title.x=element_text(size=20,face="bold"),
	       axis.title.y.right=element_text(size=20)
	      ) + xlab(TeX(xLabel)) 
	
	pt = p + annotation_custom(grob = textGrob(label = "1 block", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.01,ymax = 0.01,  xmin = -.5,xmax = 0.1) +
	annotation_custom(grob = textGrob(label = "2 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.25, ymax = 0.25,  xmin = -.5,xmax = 0.1) +
	annotation_custom(grob = textGrob(label = "3 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.65,ymax = 0.65, xmin = -.5,xmax = 0.1) +
	annotation_custom(grob = textGrob(label = "4 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.85,ymax = 0.85,  xmin = -.5,xmax = 0.1) +
	annotation_custom(grob = textGrob(label = "5 blocks", hjust = 0, gp = gpar(cex =0.9)), ymin = 0.99,ymax = 0.99,  xmin = -.5,xmax = 0.1)


	gt <- ggplot_gtable(ggplot_build(pt))
	gt$layout$clip[gt$layout$name == "panel"] <- "off"
	grid_p = grid.arrange(gt)
}
