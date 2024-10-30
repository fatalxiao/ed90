if (!require("cir")) {
  install.packages("cir")
}
library(cir)

# Interesting run (#664) from a simulated up-and-down ensemble:
# (x will be auto-generated as dose levels 1:5)
# dat=doseResponse(y=c(0,0.375,0.3846154,0.8,0.75,1),wt=c(3,8,13,10,4,1))

pieb <- read.csv("./loading volume data.csv", 1, encoding='UTF-8')
pieb

bhamou03ropi = DRtrace(x=pieb$doseSequence, y=pieb$responseSequence)
bhamou03ropi

# bhamou03ropi = DRtrace(x=c(rep(25,7),rep(22.5,2),rep(25,2),rep(27.5,6),30,rep(32.5,3),rep(35,2),rep(37.5,4),rep(40,5),rep(37.5,3),rep(40,4),42.5,rep(45,9),47.5),
#                        y=c(rep(1,8),0,1,0,rep(1,5),rep(0,2),rep(1,2),0,1,0,rep(1,3),0,rep(1,7),0,rep(1,3),rep(0,2),rep(1,8),0,1))
# bhamou03ropi

# bhamou03ropi = DRtrace(x=c(rep(25,7),rep(22.5,1),rep(25,1),rep(27.5,2),rep(30,5),rep(32.5,6),rep(30,1),rep(27.5,1),rep(30,11),rep(27.5,1),rep(30,7),rep(27.5,1),rep(30,1)),
#                        y=c(rep(1,7),rep(0,2),rep(1,1),rep(0,1),rep(1,4),rep(0,1),rep(1,7),rep(0,1),rep(1,11),rep(0,1),rep(1,7),rep(0,1),rep(1,1)))
# bhamou03ropi

# The experiment's goal is to find the 30th percentile
inv1=quickInverse(bhamou03ropi, target=0.9, adaptiveShrink = TRUE, adaptiveCurve=TRUE, conf=0.9)
inv1
# With old PAVA as the forward estimator, and without the adaptive-design corrections:
# inv0=quickInverse(dat, target=0.3, estfun=oldPAVA)


### Showing the data and the estimates
par(mar=c(3,3,1,1), mgp=c(2,.5,0), tcl=-0.25)
plot(bhamou03ropi, ylim=c(0.05,0.55), las=1) # uses plot.doseResponse()

# The true response function; true target is where it crosses the y=0.3 line
lines(seq(1,5,0.1),pweibull(seq(1,5,0.1),shape=1.1615,scale=8.4839),col=4)
abline(h=0.9,col=2,lty=3)
# Plotting the point estimates, as "tick" marks on the y=0.3 line
lines(rep(inv1$point,2),c(0.25,0.35), lwd=1.5) # CIR
# lines(rep(inv0$point,2),c(0.25,0.35),lty=2, lwd=1.5) # IR
# You could plot the CIs too,
# Here's code to plot the CIR 90\% CI as a light-green rectangle:
rect(inv1$lower90conf,0.25,inv1$upper90conf,0.35,col=rgb(0,1,0,alpha=0.3),border=NA)
#  Intervals are plotted and interval options are explored more extensively
#       in the 'deltaInverse' help page.

legend('topleft',pch=c(NA,'X',NA,NA),lty=c(1,NA,2,1),col=c(4,1,1,1),
	legend=c('True Curve','Observations','IR Estimate','CIR Estimate'),bty='n')
