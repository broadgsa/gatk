#########################################################################
# this script generates a plot of sample depth of coverage over the MHC.  
# It's rather specific to that use case, but is a good example of getting 
# Loess curve generation to work given a X/Y dataset.
# 
# 12/9/2009
# -Aaron
#########################################################################

# setup our output PNG
png(filename="bySampleJPName.png",width=1500,height=700,bg="white")

# input our data set
tbl <- read.csv("docOutJP.csv",header=TRUE) # doc_JP_SN_totalled_clean.csv

par(las=1) # make all labels horizontal
par(xpd=T, mar=par()$mar+c(0,0,-2,4)) # adjust the margins to accommodate our legend

# do the initial plot of one column of data
plot(tbl[,1],tbl[,5],xlim=c(18517983,41461957),ylim=c(0,7),type="p",cex=0.2,axes=F,ylab="Average Read Depth Of Coverage",xlab="MHC Location",col=rgb(0,0,0,0.1))

# add the custom x and y axis, so we can control their layout
axis(1,pos=0,at=seq(18517983,42061957,by=500000),col.axis="black")
axis(2,pos=18517983,at=seq(0,7,by=1),col="black")

# setup two color schemes, both with the same colors.  One has an alpha of 0.08 for the background points,
# and the other is alpha=1 for the lines (which we want to be vibrant in the foreground)
myColors <- rainbow(30,alpha=0.08)
myColors2 <- rainbow(30)

# add a legend. There is a better way to do this besides hard-coding it, but it wouldn't render correctly on my machine
legend(x=41000000,y=5,c("NA18940","NA18942","NA18943","NA18944","NA18945","NA18947","NA18948","NA18949","NA18951","NA18952","NA18953","NA18956","NA18959","NA18960","NA18961","NA18964","NA18965","NA18967","NA18968","NA18969","NA18970","NA18971","NA18972","NA18973","NA18974","NA18975","NA18976","NA18980","NA18981","NA19005"),horiz=FALSE,lty=c(1),col=c(myColors2),cex=0.8)

# loop over the remaining data sets, adding first the points to the graph, then calculating the loess points, and finally combining the points into a line
# the loess smoothing parts were inspired by: http://research.stowers-institute.org/efg/R/Statistics/loess.htm
# adjust the span value to adjust the sensitivity of curve to the local fit.
for (i in 4:33) {
    points(tbl[,1],tbl[,i],col=myColors[i],cex=0.2)
    y.loess <- loess(y ~ x, span=0.05, data.frame(x=tbl[,1], y=tbl[,i]))
    y.predict <- predict(y.loess, data.frame(x=tbl[,1]))
    lines(tbl[,1],y.predict,col=myColors2[i])
}

# close our png
dev.off()
