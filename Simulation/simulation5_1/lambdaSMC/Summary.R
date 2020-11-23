rm(list=ls())
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(deSolve)

source("ODEsimulator.R")

RMSE1_FB <- rep(NA, 40)
RMSE1_01 <- rep(NA, 40)
RMSE1_1 <- rep(NA, 40)
RMSE1_10 <- rep(NA, 40)
RMSE1_100 <- rep(NA, 40)


RMSE2_FB <- rep(NA, 40)
RMSE2_01 <- rep(NA, 40)
RMSE2_1 <- rep(NA, 40)
RMSE2_10 <- rep(NA, 40)
RMSE2_100 <- rep(NA, 40)


rep <- 40


kappa1_FB <- unlist(read.table("SMCFullB1/kappa1.txt", header = FALSE))
kappa1_01 <- unlist(read.table("SMCl01/kappa1.txt", header = FALSE))
kappa1_1 <- unlist(read.table("SMCl1/kappa1.txt", header = FALSE))
kappa1_10 <- unlist(read.table("SMCl10/kappa1.txt", header = FALSE))
kappa1_100 <- unlist(read.table("SMCl100/kappa1.txt", header = FALSE))

#kappa1_1000 <- unlist(read.table("SMCl1000/kappa1.txt", header = FALSE))
#kappa1_10000 <- unlist(read.table("SMCl10000/kappa1.txt", header = FALSE))
#kappa1 <- rbind(kappa1_FB, kappa1_01, kappa1_1, kappa1_10, kappa1_100, kappa1_1000)
value <- c(kappa1_FB, kappa1_01, kappa1_1, kappa1_10, kappa1_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
kappa1 <- data.frame(value = value, lambda = lambda)

kappa2_FB <- unlist(read.table("SMCFullB1/kappa2.txt", header = FALSE))
kappa2_01 <- unlist(read.table("SMCl01/kappa2.txt", header = FALSE))
kappa2_1 <- unlist(read.table("SMCl1/kappa2.txt", header = FALSE))
kappa2_10 <- unlist(read.table("SMCl10/kappa2.txt", header = FALSE))
kappa2_100 <- unlist(read.table("SMCl100/kappa2.txt", header = FALSE))
#kappa1_1000 <- data.frame(read.table("SMCl1000/kappa1.txt", header = FALSE))
#kappa1_10000 <- read.table("SMCl10000/kappa1.txt", header = FALSE)
#kappa2 <- rbind(kappa2_FB, kappa2_01, kappa2_1, kappa2_10, kappa2_100)
value <- c(kappa2_FB, kappa2_01, kappa2_1, kappa2_10, kappa2_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
kappa2 <- data.frame(value = value, lambda = lambda)

sigma1_FB <- unlist(read.table("SMCFullB1/sigma1.txt", header = FALSE))
sigma1_01 <- unlist(read.table("SMCl01/sigma1.txt", header = FALSE))
sigma1_1 <- unlist(read.table("SMCl1/sigma1.txt", header = FALSE))
sigma1_10 <- unlist(read.table("SMCl10/sigma1.txt", header = FALSE))
sigma1_100 <- unlist(read.table("SMCl100/sigma1.txt", header = FALSE))
#kappa1_1000 <- data.frame(read.table("SMCl1000/kappa1.txt", header = FALSE))
#kappa1_10000 <- read.table("SMCl10000/kappa1.txt", header = FALSE)
#sigma1 <- rbind(sigma1_FB, sigma1_01, sigma1_1, sigma1_10, sigma1_100)
value <- c(sigma1_FB, sigma1_01, sigma1_1, sigma1_10, sigma1_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
sigma1 <- data.frame(value = value, lambda = lambda)

sigma2_FB <- unlist(read.table("SMCFullB1/sigma2.txt", header = FALSE))
sigma2_01 <- unlist(read.table("SMCl01/sigma2.txt", header = FALSE))
sigma2_1 <- unlist(read.table("SMCl1/sigma2.txt", header = FALSE))
sigma2_10 <- unlist(read.table("SMCl10/sigma2.txt", header = FALSE))
sigma2_100 <- unlist(read.table("SMCl100/sigma2.txt", header = FALSE))
#kappa1_1000 <- data.frame(read.table("SMCl1000/kappa1.txt", header = FALSE))
#kappa1_10000 <- read.table("SMCl10000/kappa1.txt", header = FALSE)
#sigma2 <- rbind(sigma2_FB, sigma2_01, sigma2_1, sigma2_10, sigma2_100)
value <- c(sigma2_FB, sigma2_01, sigma2_1, sigma2_10, sigma2_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
sigma2 <- data.frame(value = value, lambda = lambda)

curve1_FB <- read.table("SMCFullB1/curve1.txt", header = FALSE)
curve1_01 <- read.table("SMCl01/curve1.txt", header = FALSE)
curve1_1 <- read.table("SMCl1/curve1.txt", header = FALSE)
curve1_10 <- read.table("SMCl10/curve1.txt", header = FALSE)
curve1_100 <- read.table("SMCl100/curve1.txt", header = FALSE)

for(i in 1:rep){
  RMSE1_FB[i] <- sqrt(sum((trueODE1 - curve1_FB[i,])^2)/length(trueODE1))
  RMSE1_01[i] <- sqrt(sum((trueODE1 - curve1_01[i,])^2)/length(trueODE1))
  RMSE1_1[i] <- sqrt(sum((trueODE1 - curve1_1[(i),])^2)/length(trueODE1))
  RMSE1_10[i] <- sqrt(sum((trueODE1 - curve1_10[i,])^2)/length(trueODE1))
  RMSE1_100[i] <- sqrt(sum((trueODE1 - curve1_100[i,])^2)/length(trueODE1))

}

#for(i in 1:38){
#  RMSE1_01[i] <- sqrt(sum((trueODE1 - curve1_01[i,])^2)/length(trueODE1))
#}

value <- c(RMSE1_FB, RMSE1_01, RMSE1_1, RMSE1_10, RMSE1_100)
#value <- c(curve1_FB, curve1_01, curve1_1, curve1_10, curve1_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
curve1 <- data.frame(value = value, lambda = lambda)

curve2_FB <- read.table("SMCFullB1/curve2.txt", header = FALSE)
curve2_01 <- read.table("SMCl01/curve2.txt", header = FALSE)
curve2_1 <- read.table("SMCl1/curve2.txt", header = FALSE)
curve2_10 <- read.table("SMCl10/curve2.txt", header = FALSE)
curve2_100 <- read.table("SMCl100/curve2.txt", header = FALSE)

for(i in 1:rep){
  RMSE2_FB[i] <- sqrt(sum((trueODE2 - curve2_FB[i,])^2)/length(trueODE2))
  RMSE2_01[i] <- sqrt(sum((trueODE2 - curve2_01[i,])^2)/length(trueODE2))
  RMSE2_1[i] <- sqrt(sum((trueODE2 - curve2_1[(i),])^2)/length(trueODE2))
  RMSE2_10[i] <- sqrt(sum((trueODE2 - curve2_10[i,])^2)/length(trueODE2))
  RMSE2_100[i] <- sqrt(sum((trueODE2 - curve2_100[i,])^2)/length(trueODE2))
  
}

#for(i in 1:38){
#  RMSE2_01[i] <- sqrt(sum((trueODE2 - curve2_01[i,])^2)/length(trueODE2))
#}

value <- c(RMSE2_FB, RMSE2_01, RMSE2_1, RMSE2_10, RMSE2_100)
#value <- c(curve2_FB, curve2_01, curve2_1, curve2_10, curve2_100)
lambda <- rep(c('BS', '0.1', '1', '10', '100'), times = c(40, 40, 40, 40, 40))
curve2 <- data.frame(value = value, lambda = lambda)


p0 <- ggplot(kappa1, aes(lambda, value))
p1 <- ggplot(kappa2, aes(lambda, value))
p2 <- ggplot(sigma1, aes(lambda, value))
p3 <- ggplot(sigma2, aes(lambda, value))
p4 <- ggplot(curve1, aes(lambda, value))
p5 <- ggplot(curve2, aes(lambda, value))




gname = c("lambda.eps",sep="")  
postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(theta[1]))+
            geom_hline(yintercept=2, linetype="dotted"),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(theta[2]))+
            geom_hline(yintercept=1, linetype="dotted") + ylim(c(0.8, 1.2)),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(sigma[1]))+
            geom_hline(yintercept=1, linetype="dotted"),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(sigma[2]))+
            geom_hline(yintercept=3, linetype="dotted"),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(X[1]))),
          p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = lambda))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(X[2]))),
          #         p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = Method))+rremove("x.text")+rremove("ylab")+ xlab('30')+ theme_bw()
          ncol = 6, nrow = 1, common.legend = TRUE)

dev.off()


lambda_FB <- unlist(read.table("SMCFullB1/lambda.txt", header = FALSE)[1,])
df <- data.frame(lambda = lambda_FB)
p6 <- ggplot(df, aes(x=lambda)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.02, colour="black", fill="white") + theme_bw()+rremove("ylab")+ xlab(expression(lambda))+ geom_vline(aes(xintercept=mean(lambda)), color="blue", linetype="dashed", size=1)

gname = c("lambda_sample.eps",sep="")  
postscript(gname,width=3,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(p6,
          ncol = 1, nrow = 1, common.legend = TRUE)
dev.off()


curve1_FB <- data.frame(read.table("SMCFullB/curve1.txt", header = FALSE))
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
plot(1:(length(curve1_FB[1,])),curve1_FB[1,], type = 'l', col = 1)
for(i in 2:40){
  lines(1:(length(curve1_FB[i,])),curve1_FB[i,], col = i)
}


