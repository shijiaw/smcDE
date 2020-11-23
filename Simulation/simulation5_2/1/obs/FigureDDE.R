library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

set.seed(364)
mu_m_particles = abs(rnorm(NP, 0.04, 0.02))
mu_p_particles = abs(rnorm(NP, 0.04, 0.02))
p_0_particles = rnorm(NP, 90, 20)
tau_particles = rnorm(NP, 30, 10)


sigma1 <- runif(NP, 0.3, 3)
sigma2 <- runif(NP, 0.2, 8)

gname = c("ESS.eps",sep="")  
postscript(gname,width=5,height=4,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(ESSVec, xlab = 'Time', ylab = 'ESS', type = 'b', col = 4, pch = 18)

dev.off()


gname = c("DDEparameters.eps",sep="")  
postscript(gname,width=8,height=8,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(cbind(mu_m_particles, mu_p_particles), col = 1, xlab = expression(mu[m]), ylab = expression(mu[p]), pch = 20, xlim = c(0,0.08), ylim = c(0, 0.08))
points(cbind(parametersStore[[40]][,1], parametersStore[[40]][,2]), col = 2, pch = 20)
points(cbind(parametersStore[[80]][,1], parametersStore[[80]][,2]), col = 3, pch = 20)
points(cbind(parametersStore[[120]][,1], parametersStore[[120]][,2]), col = 4, pch = 20)
points(cbind(parametersStore[[160]][,1], parametersStore[[160]][,2]), col = 5, pch = 20)
points(cbind(parametersStore[[271]][,1], parametersStore[[271]][,2]), col = 6, pch = 20)

plot(cbind(sigmaStore[[1]][,1], sigmaStore[[1]][,2]), col = 1, xlab = expression(sigma[1]), ylab = expression(sigma[2]), pch = 20, ylim = c(0,12), xlim = c(0, 18))
points(cbind(sigmaStore[[40]][,1], sigmaStore[[40]][,2]), col = 2, pch = 20)
points(cbind(sigmaStore[[80]][,1], sigmaStore[[80]][,2]), col = 3, pch = 20)
points(cbind(sigmaStore[[120]][,1], sigmaStore[[120]][,2]), col = 4, pch = 20)
points(cbind(sigmaStore[[160]][,1], sigmaStore[[160]][,2]), col = 5, pch = 20)
points(cbind(sigmaStore[[271]][,1], sigmaStore[[271]][,2]), col = 6, pch = 20)
legend("topright",c("t=1","t=40","t=80","t=120","t=160", "t=271"),cex=1.2, 
       col=1:6,pch = 20,bty="n")

plot(cbind(p_0_particles, tau_particles), col = 1, xlab = expression(p[0]), ylab = expression(tau), pch = 20, xlim = c(0,180), ylim = c(0, 80))
points(cbind(parametersStore[[40]][,3], parametersStore[[40]][,4]), col = 2, pch = 20)
points(cbind(parametersStore[[80]][,3], parametersStore[[80]][,4]), col = 3, pch = 20)
points(cbind(parametersStore[[120]][,3], parametersStore[[120]][,4]), col = 4, pch = 20)
points(cbind(parametersStore[[160]][,3], parametersStore[[160]][,4]), col = 5, pch = 20)
points(cbind(parametersStore[[271]][,3], parametersStore[[271]][,4]), col = 6, pch = 20)

dev.off()



gname = c("EstimatedDDE.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(output[,1:2], type = 'l', ylim = c(0, 25), xlab = expression(X[1]), ylab = '', lwd = 3)
lines(output[,1],fitted1, lty = 2, col = 2, lwd = 3)
lines(output[,1],upper1, lty = 3, col =2, lwd = 3)
lines(output[,1],lower1, lty = 3, col =2, lwd = 3)

plot(output[,c(1,3)], type = 'l', ylim = c(-600, 1500), xlab = expression(X[2]), ylab = '', lwd = 3)
lines(output[,1],fitted2,lty = 2, col = 2, lwd = 3)
lines(output[,1], upper2, lty = 3, col =2, lwd = 3)
lines(output[,1],lower2, lty = 3, col =2, lwd = 3)

legend("topright",c("TRUE","MEAN","CI"),cex=1.2, 
       col=c(1,2,2),lty = c(1, 2, 3), lwd = 3,pch = 20,bty="n")

dev.off()

