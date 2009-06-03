#!/broad/tools/apps/R-2.6.0/bin/Rscript

args<- commandArgs(TRUE)
verbose<- TRUE

input_fileroot <- args[1]
output_fileroot <- args[2]
dinuc <- args[3]

currentin=paste(input_fileroot,dinuc,"csv",sep=".");
currentout=paste(output_fileroot,dinuc,"parameters",sep=".");
con <- file(currentin, "r")
con_out <- file(currentout, "w+")
data<-read.table(con,header=TRUE,sep=",")
reg<-glm(indicator~1+logitQ+I(logitQ^2)+I(logitQ^3)+I(logitQ^4)+pos+I(logitQ*pos)+I((logitQ^2)*pos)+I((logitQ^3)*pos)+I((logitQ^4)*pos)+I(pos^2)+I(logitQ*(pos^2))+I((logitQ^2)*(pos^2))+I((logitQ^3)*(pos^2))+I((logitQ^4)*(pos^2))+I(pos^3)+I(logitQ*(pos^3))+I((logitQ^2)*(pos^3))+I((logitQ^3)*(pos^3))+I((logitQ^4)*(pos^3))+I(pos^4)+I(logitQ*(pos^4))+I((logitQ^2)*(pos^4))+I((logitQ^3)*(pos^4))+I((logitQ^4)*(pos^4)),family=binomial("logit"),data=data, weights=data$count)
write(coefficients(reg),file=con_out)
