###Plot of ODE and DDE ###
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)

parameter1 <- c(thetaStore[[1]][,1], thetaStore[[203]][,1], thetaStore[[610]][,1], thetaStore[[1220]][,1])
 
parameter2 <- c(thetaStore[[1]][,2], thetaStore[[203]][,2], thetaStore[[610]][,2], thetaStore[[1220]][,2])

data1 <- data.frame(kappa1 = parameter1, kappa2 = parameter2, iter = rep(c("1","2","3","4"), each = 500))


p1 <- ggplot(data1, aes(x=kappa1, y=kappa2, color=iter)) +
  geom_point()+ scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "blue"))+xlim(c(-5, 5))+ylim(c(-1,5)) +
  theme_bw()+ xlab(expression(theta[1]))+ ylab(expression(theta[2]))+
  geom_point(aes(x=-2, y=1), colour="black")+geom_point(aes(x=2, y=1), colour="black")


gname = c("SMCnoBasis.eps",sep="")  
postscript(gname,width=4,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
  p1,
  ncol = 1, nrow = 1, labels = NULL,legend = 'none'
)
dev.off()
