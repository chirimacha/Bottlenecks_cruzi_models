makeAllPlots <- function(doPrev=TRUE, onSameAxes=FALSE, letterLabels=FALSE,
                         showRepresentative=TRUE) {

    # Do plots for all bottleneck lengths: 0, 14, 30, 60, 90, 180 days
    lens <- c(0, 14, 30, 60, 90, 180)
    bStart <- 120
    #highRunIndex <- c(25, 206, 136, 67, 54, 65)
    #medRunIndex <- c(11, 23, 137, 37, 61, 69)
    #lowRunIndex <- c(16, 10, 138, 38, 58, 59)
    labels <- c("A", "B", "C", "D", "E", "F")

    if (onSameAxes) {
        #X11(display="", 10, 7)
        par(mfrow=c(2, 3))
    }
    
    for (i in 1:length(lens)) {
        fName <- paste('~/Dropbox/Superspreading/New Bug-Cuy Code/aaron_output/bottleneck_length_1000_reps/vector_prev_vs_time_bottleLength_',lens[i],'.csv',sep='')
        vectorPrev <<- read.csv(fName)[-1]
        if (letterLabels)
            labText <- labels[i]
        else
            labText <- ""
        makePlot(doPrev, lens[i], bStart, label=labText, showRepresentative)
    }
    
    par(mfrow=c(1,1))
    
    # Add the x and y labels:
    title(xlab="days since start of simulation", line=3.5, cex.lab=1.2)
    if (doPrev)
        title(ylab="infection prevalence in vectors (%)", line=3, cex.lab=1.2)
    else
        title(ylab="number of infected vectors", line=3, cex.lab=1.2)

}

# highRun, medRun, and lowRun are the indices of runs that have high,
# medium, and low prevalence.
makePlot <- function(doPrev = FALSE, lengthBottle, startBottle,
                     label="", showRepresentative=TRUE) {
    
    if (!exists("vectorPrev"))
        vectorPrev <<- read.csv("~/Dropbox/Superspreading/New Bug-Cuy Code/output/cuybottle_nobugbottle_1000reps/vector_prev_vs_time.csv")[-1]

    simLengthDays <- ncol(vectorPrev)
    
    # Thin the data: show only 10% of the runs (100 total runs)
    x <- 1:1000
    runs <- sample(x, 100)
    runData <- vectorPrev[runs,]
    nsims <- 100
    
    # default is to show the number of bugs infected on y axis:
    nbugs <- 1000
    ylimit <- nbugs
    divFactor <- 1
    if (doPrev) {
        # change the defaults to show prevalence (%) on y axis:
        ylimit <- 100
        divFactor <- 10
    }

    #color<-rainbow(nsims)
    color <- gray.colors(nsims)
    
    plot(1:simLengthDays, runData[1,]/divFactor, xlab="", ylab="", ylim=c(0, ylimit), type="l", cex.axis=1.2)
    mtext(label, side=3, adj=0, line=1.2, font=2)

    for(phi in 2:nsims)
        points(1:simLengthDays, runData[phi,]/divFactor, type="l", col=color[phi])
        
    # Show the bottleneck start and end times with vertical lines
    if (lengthBottle > 0) {
        shifts<-seq(0,simLengthDays,365)
        abline(v=c(startBottle+shifts),lty=1, lwd=2)
        abline(v=c(startBottle+lengthBottle+shifts),lty=2, lwd=2)
    }

    if (showRepresentative) {
        # Add a few representative runs in color. Use the runs with the
        # 5th, 50th, and 95th percentile of area under the curve
        # (i.e., total sum of vectors infected).
        sumInfected <- rowSums(runData)
        quants <- quantile(sumInfected, c(.05, .5, .95))
        colorToPlot <- c("darkgreen", "blue", "red")
        for (j in 1:3) {
            cat("doing j=",j,"\n")
            absTest <- abs(sumInfected - quants[j])
            ind <- which(absTest == min(absTest))
            cat("ind=",ind,"\n")
            points(1:simLengthDays, runData[ind[1],]/divFactor, type="l", col=colorToPlot[j], lwd=2)
        }
    }
}

getIndices <- function() {
    for (i in 1:1000) {
        if (max(vectorPrev[i,]) > 800)
            cat(i,'\n')
    }
}
