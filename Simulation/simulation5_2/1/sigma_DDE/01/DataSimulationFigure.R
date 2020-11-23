###Plot of ODE and DDE ###
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

ODEoutput <- data.frame(
  time = output[,1],
  X1 = output[,2],
  Y1 = y1
)

gname = c("DDEtrajectory.eps",sep="")  
postscript(gname,width=4,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
ggplot(ODEoutput, aes(time)) + 
  geom_line(aes(y = X1, colour = "var0")) + 
  geom_point(aes(y = Y1, colour = "var1"))+ theme_bw()+rremove("ylab"),
ncol = 1, nrow = 1, labels = NULL,legend = 'none'
)
dev.off()


