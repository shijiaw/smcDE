rm(list=ls())
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(deSolve)

source("ODEsimulator.R")


RMSE1_1 <- rep(NA, 40)
RMSE1_2 <- rep(NA, 40)
RMSE1_4 <- rep(NA, 40)
RMSE1_6 <- rep(NA, 40)
RMSE1_10 <- rep(NA, 40)

RMSE2_1 <- rep(NA, 40)
RMSE2_2 <- rep(NA, 40)
RMSE2_4 <- rep(NA, 40)
RMSE2_6 <- rep(NA, 40)
RMSE2_10 <- rep(NA, 40)


rep <- 40

kappa1_1 <- unlist(read.table("SMC1/kappa1.txt", header = FALSE))
kappa1_2 <- unlist(read.table("SMC2/kappa1.txt", header = FALSE))
kappa1_4 <- unlist(read.table("SMC4/kappa1.txt", header = FALSE))
kappa1_6 <- unlist(read.table("SMC6/kappa1.txt", header = FALSE))
kappa1_10 <- unique(unlist(read.table("SMC10/kappa1.txt", header = FALSE)))

value <- c(kappa1_10, kappa1_6, kappa1_4, kappa1_2, kappa1_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

kappa1 <- data.frame(value = value, nknots = nknots)

kappa2_1 <- unlist(read.table("SMC1/kappa2.txt", header = FALSE))
kappa2_2 <- unlist(read.table("SMC2/kappa2.txt", header = FALSE))
kappa2_4 <- unlist(read.table("SMC4/kappa2.txt", header = FALSE))
kappa2_6 <- unlist(read.table("SMC6/kappa2.txt", header = FALSE))
kappa2_10 <- unique(unlist(read.table("SMC10/kappa2.txt", header = FALSE)))

value <- c(kappa2_10, kappa2_6, kappa2_4, kappa2_2, kappa2_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

kappa2 <- data.frame(value = value, nknots = nknots)


sigma1_1 <- unlist(read.table("SMC1/sigma1.txt", header = FALSE))
sigma1_2 <- unlist(read.table("SMC2/sigma1.txt", header = FALSE))
sigma1_4 <- unlist(read.table("SMC4/sigma1.txt", header = FALSE))
sigma1_6 <- unlist(read.table("SMC6/sigma1.txt", header = FALSE))
sigma1_10 <- unlist(read.table("SMC10/sigma1.txt", header = FALSE))

value <- c(sigma1_10, sigma1_6, sigma1_4, sigma1_2, sigma1_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

sigma1 <- data.frame(value = value, nknots = nknots)

sigma2_1 <- unlist(read.table("SMC1/sigma2.txt", header = FALSE))
sigma2_2 <- unlist(read.table("SMC2/sigma2.txt", header = FALSE))
sigma2_4 <- unlist(read.table("SMC4/sigma2.txt", header = FALSE))
sigma2_6 <- unlist(read.table("SMC6/sigma2.txt", header = FALSE))
sigma2_10 <- unique(unlist(read.table("SMC10/sigma2.txt", header = FALSE)))

value <- c(sigma2_10, sigma2_6, sigma2_4, sigma2_2, sigma2_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

sigma2 <- data.frame(value = value, nknots = nknots)

curve1_1 <- read.table("SMC1/curve1.txt", header = FALSE)
curve1_2 <- read.table("SMC2/curve1.txt", header = FALSE)
curve1_4 <- read.table("SMC4/curve1.txt", header = FALSE)
curve1_6 <- read.table("SMC6/curve1.txt", header = FALSE)
curve1_10 <- read.table("SMC10/curve1.txt", header = FALSE)


for(i in 1:rep){
  RMSE1_1[i] <- sqrt(sum((trueODE1 - curve1_1[i,])^2)/length(trueODE1))
  RMSE1_2[i] <- sqrt(sum((trueODE1 - curve1_2[i,])^2)/length(trueODE1))
  RMSE1_4[i] <- sqrt(sum((trueODE1 - curve1_4[i,])^2)/length(trueODE1))
  RMSE1_6[i] <- sqrt(sum((trueODE1 - curve1_6[i,])^2)/length(trueODE1))
  RMSE1_10[i] <- sqrt(sum((trueODE1 - curve1_10[i,])^2)/length(trueODE1))
}

value <- c(RMSE1_10, RMSE1_6, RMSE1_4, RMSE1_2, RMSE1_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

curve1 <- data.frame(value = value, nknots = nknots)


curve2_1 <- read.table("SMC1/curve2.txt", header = FALSE)
curve2_2 <- read.table("SMC2/curve2.txt", header = FALSE)
curve2_4 <- read.table("SMC4/curve2.txt", header = FALSE)
curve2_6 <- read.table("SMC6/curve2.txt", header = FALSE)
curve2_10 <- read.table("SMC10/curve2.txt", header = FALSE)


for(i in 1:rep){
  RMSE2_1[i] <- sqrt(sum((trueODE2 - curve2_1[i,])^2)/length(trueODE1))
  RMSE2_2[i] <- sqrt(sum((trueODE2 - curve2_2[i,])^2)/length(trueODE1))
  RMSE2_4[i] <- sqrt(sum((trueODE2 - curve2_4[i,])^2)/length(trueODE1))
  RMSE2_6[i] <- sqrt(sum((trueODE2 - curve2_6[i,])^2)/length(trueODE1))
  RMSE2_10[i] <- sqrt(sum((trueODE2 - curve2_10[i,])^2)/length(trueODE1))
}

value <- c(RMSE2_10, RMSE2_6, RMSE2_4, RMSE2_2, RMSE2_1)
nknots <- rep(c('7', '11', '16', '31', '61'), each = 40)
nknots <- as.factor(nknots)
nknots <- factor(nknots,levels(nknots)[c(5,1,2,3,4)])

curve2 <- data.frame(value = value, nknots = nknots)



p0 <- ggplot(kappa1, aes(nknots, value))
p1 <- ggplot(kappa2, aes(nknots, value))
p2 <- ggplot(sigma1, aes(nknots, value))
p3 <- ggplot(sigma2, aes(nknots, value))
p4 <- ggplot(curve1, aes(nknots, value))
p5 <- ggplot(curve2, aes(nknots, value))




gname = c("nknots.eps",sep="")  
postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)

ggarrange(p0 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(theta[1]))+
            geom_hline(yintercept=2, linetype="dotted"),
          p1 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(theta[2]))+
            geom_hline(yintercept=1, linetype="dotted") + ylim(c(0.8, 1.2)),
          p2 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(sigma[1]))+
            geom_hline(yintercept=1, linetype="dotted"),
          p3 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(sigma[2]))+
            geom_hline(yintercept=3, linetype="dotted"),
          p4 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(X[1]))),
          p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = nknots))+ theme_bw()+rremove("x.text")+rremove("ylab")+ xlab(expression(RMSE(X[2]))),
          #         p5 + geom_boxplot(fill = "white", colour = "#3366FF", outlier.colour = "red", outlier.shape = 1) + geom_boxplot(aes(color = Method))+rremove("x.text")+rremove("ylab")+ xlab('30')+ theme_bw()
          ncol = 6, nrow = 1, common.legend = TRUE)

dev.off()




