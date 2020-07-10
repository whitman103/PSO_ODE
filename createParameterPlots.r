library(extrafont)
library(Cairo)
graphics.off()
CairoPDF("outPlot.pdf",10,7.5)
plot.new()
trueData=read.table("threeFolder//trueOutPlot.txt",header=FALSE)
fitData=read.table("threeFolder//fitData.txt",header=FALSE)
par(family="CMU Serif")
numDecay=5
numPlots=(length(trueData)-numDecay)/2+numDecay
numParameterPlots=12
options(digits=3)
par(mfrow=c(3,numParameterPlots/3),mai=c(.6,.5,.15,.2),omi=c(.1,.1,.1,.1))
for(plot in 0:(numParameterPlots-1)){
    #xMin=min(c(trueData[,2*plot+1],fitData[,2*plot+1]))
    #xMax=max(c(trueData[,2*plot+1],fitData[,2*plot+1]))
    #yMin=min(c(trueData[,2*plot+2],fitData[,2*plot+2]))
    #yMax=max(c(trueData[,2*plot+2],fitData[,2*plot+2]))
    xMin=0
    xMax=50
    yMin=0
    yMax=0.05
    xRange=xMax-xMin
    xMin=xMin-.1*xRange
    xMax=xMax+.1*xRange
    yRange=yMax-yMin
    yMin=yMin-.1*yRange
    yMax=yMax+.1*yRange
    fitPower=round(mean(fitData[,3*plot+3]),3)
    truePower=round(trueData[,3*plot+3],3)
    plot(fitData[,3*plot+1],fitData[,3*plot+2],xlim=c(xMin,xMax),ylim=c(yMin,yMax),xlab=paste("Const",plot),ylab=paste("Norm",plot),pch=19,cex=.75,main=paste("fitP=",fitPower," trueP=",truePower))
    points(trueData[,3*plot+1],trueData[,3*plot+2],col="red",pch=19,cex=2)
    points(mean(fitData[,3*plot+1]),mean(fitData[,3*plot+2]),col="blue",pch=19,cex=2)
}
dev.off()
