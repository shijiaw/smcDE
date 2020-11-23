library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

set.seed(364)
kappa1 <- runif(NP, 0, 4)
kappa2 <- runif(NP, 0, 4)

sigma1 <- runif(NP, 0.3, 3)
sigma2 <- runif(NP, 0.3, 3)

gname = c("ESS.eps",sep="")  
postscript(gname,width=5,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(ESSVec, xlab = 'Time', ylab = 'ESS', type = 'b', col = 4, pch = 18)

dev.off()


gname = c("ParameterBasis.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(cbind(kappaStore[[1]][,1], kappaStore[[1]][,2]), col = 1, xlab = expression(kappa[1]), ylab = expression(kappa[2]), pch = 20, xlim = c(-5,5), ylim = c(-0.5, 4))
points(cbind(kappaStore[[40]][,1], kappaStore[[40]][,2]), col = 2, pch = 20)
points(cbind(kappaStore[[80]][,1], kappaStore[[80]][,2]), col = 3, pch = 20)
points(cbind(kappaStore[[120]][,1], kappaStore[[120]][,2]), col = 4, pch = 20)
points(cbind(kappaStore[[160]][,1], kappaStore[[160]][,2]), col = 5, pch = 20)
points(cbind(kappaStore[[193]][,1], kappaStore[[193]][,2]), col = 6, pch = 20)

plot(cbind(sigmaStore[[1]][,1], sigmaStore[[1]][,2]), col = 1, xlab = expression(sigma[1]), ylab = expression(sigma[2]), pch = 20, ylim = c(0,12), xlim = c(0, 18))
points(cbind(sigmaStore[[40]][,1], sigmaStore[[40]][,2]), col = 2, pch = 20)
points(cbind(sigmaStore[[80]][,1], sigmaStore[[80]][,2]), col = 3, pch = 20)
points(cbind(sigmaStore[[120]][,1], sigmaStore[[120]][,2]), col = 4, pch = 20)
points(cbind(sigmaStore[[160]][,1], sigmaStore[[160]][,2]), col = 5, pch = 20)
points(cbind(sigmaStore[[193]][,1], sigmaStore[[193]][,2]), col = 6, pch = 20)
legend("topright",c("t=1","t=40","t=80","t=120","t=160", "t=193"),cex=1.2, 
       col=1:6,pch = 20,bty="n")

dev.off()



gname = c("EstimatedODEwithLambda.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(output[-1,2], type = 'l', ylim = c(-8, 11), xlab = expression(X[1]), ylab = '', lwd = 3)
lines(fitted1, lty = 2, col = 2, lwd = 3)
lines(upper1, lty = 3, col =2, lwd = 3)
lines(lower1, lty = 3, col =2, lwd = 3)

plot(output[-1,3], type = 'l', ylim = c(-80, 100), xlab = expression(X[2]), ylab = '', lwd = 3)
lines(fitted2,lty = 2, col = 2, lwd = 3)
lines(upper2, lty = 3, col =2, lwd = 3)
lines(lower2, lty = 3, col =2, lwd = 3)

legend("topright",c("TRUE","MEAN","CI"),cex=1.2, 
       col=c(1,2,2),lty = c(1, 2, 3), lwd = 3,pch = 20,bty="n")

dev.off()

