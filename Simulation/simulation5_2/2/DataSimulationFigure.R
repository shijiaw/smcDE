###Plot of ODE and DDE ###
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

ODEoutput <- data.frame(
  time = output[,1],
  X1 = output[,2],
  X2 = output[,3],
  Y1 = y1,
  Y2 = y2
)

gname = c("DDEsimulator.eps",sep="")  
postscript(gname,width=8,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
ggplot(ODEoutput, aes(time)) + 
  geom_line(aes(y = X1, colour = "var0")) + 
  geom_point(aes(y = Y1, colour = "var1"))+ theme_bw()+rremove("ylab")+ xlab(expression(X[1](t))),
ggplot(ODEoutput, aes(time)) + 
  geom_line(aes(y = X2, colour = "var0")) + 
  geom_point(aes(y = Y2, colour = "var1"))+ theme_bw()+rremove("ylab")+ xlab(expression(X[2](t))),
ncol = 2, nrow = 1, labels = NULL,legend = 'none'
)
dev.off()


DDEoutput <- data.frame(
  time = output[,1],
  X1 = output[,2],
  X2 = output[,3],
  Y1 = y1,
  Y2 = y2
)

gname = c("DDEsimulator.eps",sep="")  
postscript(gname,width=8,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
  ggplot(DDEoutput, aes(time)) + 
    geom_line(aes(y = X1, colour = "var0")) + 
    geom_point(aes(y = Y1, colour = "var1"))+ theme_bw()+rremove("ylab")+ xlab(expression(X[1](t))),
  ggplot(DDEoutput, aes(time)) + 
    geom_line(aes(y = X2, colour = "var0")) + 
    geom_point(aes(y = Y2, colour = "var1"))+ theme_bw()+rremove("ylab")+ xlab(expression(X[2](t))),
  ncol = 2, nrow = 1, labels = NULL,legend = 'none'
)
dev.off()
