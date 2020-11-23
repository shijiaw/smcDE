gname = c("EstimatedDDE2.eps",sep="")  
postscript(gname,width=10,height=5,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,2),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

plot(output[,2], type = 'l', ylim = c(-8, 30), xlab = expression(X[1]), ylab = '', lwd = 3)
lines(fitted1, lty = 2, col = 2, lwd = 3)
lines(upper1, lty = 3, col =2, lwd = 3)
lines(lower1, lty = 3, col =2, lwd = 3)

plot(output[,3], type = 'l', ylim = c(-60, 600), xlab = expression(X[2]), ylab = '', lwd = 3)
lines(fitted2,lty = 2, col = 2, lwd = 3)
lines(upper2, lty = 3, col =2, lwd = 3)
lines(lower2, lty = 3, col =2, lwd = 3)

legend("topright",c("TRUE","MEAN","CI"),cex=1.2, 
       col=c(1,2,2),lty = c(1, 2, 3), lwd = 3,pch = 20,bty="n")

dev.off()

