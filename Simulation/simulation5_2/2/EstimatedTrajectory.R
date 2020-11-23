###Plot of ODE and DDE ###
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


curve_matrix1 <- matrix(NA, nr = nrow(basismat), nc = NP)
for(i in 1:NP){
  curve_matrix1[,i] <- basismat%*%c1_estimate[i,]
}

upper_curve1 <- rep(NA, nrow(basismat))
mean_curve1 <- rep(NA, nrow(basismat))
lower_curve1 <- rep(NA, nrow(basismat))
for(j in 1:nrow(basismat)){
  upper_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.975, type='quantile', weights=W, normwt=TRUE)
  lower_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.025, type='quantile', weights=W, normwt=TRUE)
  mean_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.5, type='quantile', weights=W, normwt=TRUE)
}


plot(mean_curve1, type = 'l', ylim = c(-8, 30))
lines(output[,2], type = 'l', col = 3)
lines(upper_curve, type = 'l', col = 2)
lines(lower_curve, type = 'l', col = 2)




curve_matrix2 <- matrix(NA, nr = nrow(basismat), nc = NP)
for(i in 1:NP){
  curve_matrix2[,i] <- basismat%*%c2_estimate[i,]
}

upper_curve2 <- rep(NA, nrow(basismat))
mean_curve2 <- rep(NA, nrow(basismat))
lower_curve2 <- rep(NA, nrow(basismat))
for(j in 1:nrow(basismat)){
  upper_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.975, type='quantile', weights=W, normwt=TRUE)
  lower_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.025, type='quantile', weights=W, normwt=TRUE)
  mean_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.5, type='quantile', weights=W, normwt=TRUE)
}

plot(mean_curve2, type = 'l', ylim = c(-60, 600))
lines(output[,3], type = 'l', col = 3)
lines(upper_curve2, type = 'l', col = 2)
lines(lower_curve2, type = 'l', col = 2)



variable1 = rep(c("X1", "50%", "97.5%", "2.5%"), each = 101)
variable2 = rep(c("X2", "50%", "97.5%", "2.5%"), each = 101)
variable1 <- as.factor(variable1)
#variable1 <- factor(variable1,levels(variable1)[c(4,1,2,3)])
variable2 <- as.factor(variable2)
#variable2 <- factor(variable2,levels(variable2)[c(4,1,2,3)])


ODEoutput1 <- data.frame(
  time = rep(output[,1], 4),
  value = c(output[,2], mean_curve1, upper_curve1, lower_curve1),
  variable = variable1
)

ODEoutput2 <- data.frame(
  time = rep(output[,1], 4),
  value = c(output[,3], mean_curve2, upper_curve2, lower_curve2),
  variable = variable2
)



p1 <- ggplot(ODEoutput1, aes(x=time, y=value, group=variable)) + geom_line(aes(linetype=variable, color=variable))+ 
  scale_linetype_manual(values=c("dotted","twodash", "dotted", "solid")) +
  scale_color_manual(values=c('blue','red','blue','black')) + theme_bw()+rremove("x.text")+rremove("ylab")+ xlab('X1') 
p2 <- ggplot(ODEoutput2, aes(x=time, y=value, group=variable)) + geom_line(aes(linetype=variable, color=variable))+ 
  scale_linetype_manual(values=c("dotted","twodash", "dotted", "solid")) +
  scale_color_manual(values=c('blue','red','blue','black')) + theme_bw()+rremove("x.text")+rremove("ylab")+ xlab('X2') 


gname = c("DDEtrajectory.eps",sep="")  
postscript(gname,width=8,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(
 p1,
 p2,
  ncol = 2, nrow = 1, labels = NULL,legend = 'none'
)
dev.off()

