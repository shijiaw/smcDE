library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)


sigma1 <- data.frame(error = rep(0.1, 12))
sigma1$W = unlist(read.table("01/W.txt", header = FALSE))
sigma1$nu = unlist(read.table("01/nu.txt", header = FALSE))
sigma1$P = unlist(read.table("01/P.txt", header = FALSE))
sigma1$tau = unlist(read.table("01/tau.txt", header = FALSE))
#sigma1$sigma = unlist(read.table("01/sigma.txt", header = FALSE))
sigma1$RMSE = unlist(read.table("01/MSE.txt", header = FALSE))

sigma5 <- data.frame(error = rep(0.5, 12))
sigma5$W = unlist(read.table("05/W.txt", header = FALSE))
sigma5$nu = unlist(read.table("05/nu.txt", header = FALSE))
sigma5$P = unlist(read.table("05/P.txt", header = FALSE))
sigma5$tau = unlist(read.table("05/tau.txt", header = FALSE))
#sigma5$sigma = unlist(read.table("05/sigma.txt", header = FALSE))
sigma5$RMSE = unlist(read.table("05/MSE.txt", header = FALSE))

sigma15 <- data.frame(error = rep(1.5, 12))
sigma15$W = unlist(read.table("15/W.txt", header = FALSE))
sigma15$nu = unlist(read.table("15/nu.txt", header = FALSE))
sigma15$P = unlist(read.table("15/P.txt", header = FALSE))
sigma15$tau = unlist(read.table("15/tau.txt", header = FALSE))
#sigma15$sigma = unlist(read.table("15/sigma.txt", header = FALSE))
sigma15$RMSE = unlist(read.table("15/MSE.txt", header = FALSE))

Results <- rbind(sigma1, sigma5, sigma15)
Results$sigma = as.factor(Results$error)

#Resultssummary <- summarySE(Results, measurevar=c("nu"), groupvars=c("sigma"))
p0 <- ggplot(Results, aes(sigma, nu))
p1 <- ggplot(Results, aes(sigma, W))
p2 <- ggplot(Results, aes(sigma, P))
p3 <- ggplot(Results, aes(sigma, tau))
p4 <- ggplot(Results, aes(sigma, RMSE))

gname = c("DDE.eps",sep="")  
postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(nu))+
            geom_hline(yintercept=0.8, linetype="dotted"),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(P))+
            geom_hline(yintercept=2, linetype="dotted"),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(tau))+
            geom_hline(yintercept=3, linetype="dotted"),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(W))+
            geom_hline(yintercept=8.16, linetype="dotted"),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = sigma))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab("RMSE"),
          ncol = 5, nrow = 1, common.legend = TRUE)

dev.off()


