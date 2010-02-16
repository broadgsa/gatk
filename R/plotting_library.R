#!/broad/tools/apps/R-2.6.0/bin/Rscript

args <- commandArgs(TRUE)
fileToRead <- args[1]
functionToRun <- args[2]
functionSpecificArgs <- args[3]

## load the function to run

if ( funtionToRun == "PlotInterleavedRows" ) {
### PLOT INTERLEAVED ROWS FUNCTION ###
#  - expects a file of the form
#
# sample_a \t 0.8 \t 0.6 \t 0.5
# sample_a \t 0 \t 1 \t 3
# sample_b \t 0.5 \t 0.3 \t 0.1
# sample_b \t 1 \t 2 \t 4
# 
# and an argument string
#   x_label;y_label;plot_title;base_name_for_pdf
#  - end of info -
### PLOT INTERLEAVED ROWS FUNCTION ###
PlotInterleavedRows <- function(inFile,args) {
	arglist = unlist(strsplit(args,";"))
	xlabel = arglist[1]
	ylabel = arglist[2]
	title = arglist[3]
	outFileBase = arglist[4]
	
	allPoints <- as.matrix(read.table(inFile))
	# set up colors
	colors = rainbow(ncol(allPoints)-1,s=0.8,v=0.8,gamma=0.6,start=0.0,end=0.9)
	styles = c(rep(1,ncol(allPoints)-1))
	evalPoints = matrix(nrow=nrow(allPoints)/2,ncol=ncol(allPoints))
	funcVal = matrix(nrow=nrow(allPoints)/2,ncol=ncol(allPoints))
	# convert to two matrices by de-interleaving and transposing
	for ( i in 1:(nrow(allPoints)/2) ) {
		evalPoints[i,] <- allPoints[2*i,]
		funcVal[i,] <- allPoints[2*i-1,]
	}
	
	evalPoints <- t(evalPoints)
	funcVal <- t(funcVal)
	# plot and put legend on
	pdf(paste(outFileBase,"_rplot",".pdf",sep=""))
	matplot(evalPoints,funcVal,col=colors,lty=styles,"l",xlab=xlabel,ylab=ylabel)
	legend("topright",funcVal[1,],lty=styles,col=colors)
	title(main=title,outer=TRUE)
	# save
	dev.off()
}

PlotInterleavedRows(fileToRead,functionSpecificArgs)

}

if ( functionToRun == "PlotHeatmap" ) {
### PLOT HEATMAP FUNCTION ###
#
# Normally what is meant by "heatmap" is just an image() of the
# matrix; in accordance with that, THIS FUNCTION DOES NOT COMPUTE
# DENDROGRAMS THROUGH HEATMAP(), so no rows and columns are not
# re-ordered, and dendrograms are not displayed.
#
# - expects a file of the form
#
# rentry1 \t rentry2 \t rentry3 \t ...
# colentry1          \t  0.7    \t   0.9   \t   0.4   \t ...
# colentry2          \t  0.8    \t   0.7   \t   0.6   \t ...
#   ...
# Note that the rows and columns don't line up. R understands this
# and deals with it.
# Also expects an argument string:
#   row_label;column_label;data_rescale_factor;plot_title;base_name_for_pdf
# - end of info -
### PLOT HEATMAP FUNCTION ###
PlotHeatmap <- function(inFile,args) {
	arglist = unlist(strsplit(args,split=";"))
	row_label = arglist[1]
	column_label = arglist[2]
	data_rescale_factor <- as.numeric(arglist[3])
	plot_title = arglist[4]
	base_name_for_pdf = arglist[5]
	image_matrix <- data_rescale_factor*as.matrix(read.table(inFile))
	## change default colors to include "cool" colors for lower end of spectrum
	## e.g. red ~ near 1, yellow ~ near .75, green ~ near .5, teal ~ near .25
	##      blue ~ near 0
	colors <- rev(rainbow(32,start=0,end=0.6,s=0.9,v=0.9,gamma=0.8))
	pdf(paste(base_name_for_pdf,"_rplot",".pdf",sep=""))
	heatmap(image_matrix,Rowv=NA,Colv=NA,ylab=row_label,xlab=column_label,col=colors)
	title(main=plot_title,outer=TRUE)
	dev.off()
}

PlotHeatmap(fileToRead,functionSpecificArgs)

}