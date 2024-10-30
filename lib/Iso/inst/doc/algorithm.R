### R code from vignette source 'algorithm.Rnw'

###################################################
### code chunk number 1: algorithm.Rnw:27-28
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: algorithm.Rnw:33-39
###################################################
library(Iso)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "Iso"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "Iso"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: algorithm.Rnw:475-483
###################################################
getOption("SweaveHooks")[["fig"]]()
require(Iso)
OP <- par(mfrow=c(3,2),mar=c(4,4,3,1))
for(i in 2:6) {
   plot(ufit(vigour[,i],x=vigour[,1]),type="l",ylim=c(0,0.3),
        xlab="year",ylab="vigour",main=paste("stand",i-1),cex.main=1.5)
   points(vigour[,1],vigour[,i],pch="+",col="red")
}
par(OP)


###################################################
### code chunk number 4: algorithm.Rnw:501-507 (eval = FALSE)
###################################################
## par(mfrow=c(3,2),mar=c(4,4,3,1))
## for(i in 2:6) {
##    plot(ufit(vigour[,i],x=vigour[,1]),type="l",ylim=c(0,0.3),
##         xlab="year",ylab="vigour",main=paste("stand",i-1),cex.main=1.5)
##    points(vigour[,1],vigour[,i],pch="+",col="red")
## }


###################################################
### code chunk number 5: algorithm.Rnw:520-526 (eval = FALSE)
###################################################
##    xm <- apply(vigour[,2:6],1,mean)
##    par(mar=c(4,4,3,1))
##    plot(ufit(xm,x=vigour[,1]),type="l",ylim=c(0,0.3),
##         xlab="year",ylab="vigour",main="Mean over stands",cex.main=1.5)
##    points(vigour[,1],xm,pch=22,col="red")
##    for(i in 2:6) points(vigour[,1],vigour[,i],pch="+",col="blue")


###################################################
### code chunk number 6: algorithm.Rnw:531-537
###################################################
getOption("SweaveHooks")[["fig"]]()
   xm <- apply(vigour[,2:6],1,mean)
   par(mar=c(4,4,3,1))
   plot(ufit(xm,x=vigour[,1]),type="l",ylim=c(0,0.3),
        xlab="year",ylab="vigour",main="Mean over stands",cex.main=1.5)
   points(vigour[,1],xm,pch=22,col="red")
   for(i in 2:6) points(vigour[,1],vigour[,i],pch="+",col="blue")


###################################################
### code chunk number 7: algorithm.Rnw:553-554
###################################################
tools::compactPDF(".",gs_quality="ebook")


