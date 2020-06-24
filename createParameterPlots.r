trueData=read.table("trueOutPlot.txt",header=FALSE)
fitData=read.table("fitData.txt",header=FALSE)

numDecay=5
numPlots=(length(trueData)-numDecay)/2+numDecay
numParameterPlots=(length(trueData)-numDecay)/2
par(mfrow=c(3,numParameterPlots/3))
for(plot in 0:(numParameterPlots-1)){
    xMin=min(c(trueData[,2*plot+1],fitData[,2*plot+1]))
    xMax=max(c(trueData[,2*plot+1],fitData[,2*plot+1]))
    yMin=min(c(trueData[,2*plot+2],fitData[,2*plot+2]))
    yMax=max(c(trueData[,2*plot+2],fitData[,2*plot+2]))
    xRange=xMax-xMin
    xMin=xMin-.1*xRange
    xMax=xMax+.1*xRange
    yRange=yMax-yMin
    yMin=yMin-.1*yRange
    yMax=yMax+.1*yRange
    plot(fitData[,2*plot+1],fitData[,2*plot+2],xlim=c(xMin,xMax),ylim=c(yMin,yMax),xlab=paste("Const",plot),ylab=paste("Norm",plot))
    points(trueData[,2*plot+1],trueData[,2*plot+2],col="red")
    points(mean(fitData[,2*plot+1]),mean(fitData[,2*plot+2]),col="blue")
}
