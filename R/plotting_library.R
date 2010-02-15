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
	# save
	dev.off()
}

PlotInterleavedRows(fileToRead,functionSpecificArgs)

}