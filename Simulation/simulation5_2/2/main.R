rm(list=ls())
library(deSolve)
library("truncnorm")
#library("fda")
library(Hmisc)
#library("MCMCpack")
library(numDeriv)
library(fda)
library(NLRoot)
library(smcUtils)
library(MASS)

NITER <- 1
knotsNpoints <- 20
#lambda <- 10
NP <- 300
lambda <- rgamma(NP, 1, 1)
CESSthresholds <- 0.3
resampleThreshold <- 0.5
seeds <- (1:NITER)
filename_mu_m = './mu_m'
filename_mu_p = './mu_p'
filename_p_0 = './p_0'
filename_tau = './tau'
filename_sigma1 = './sigma1'
filename_sigma2 = './sigma2'
filename_curve1 = './curve1'
filename_curve2 = './curve2'
filename_lambda = './lambda'

filename_mu_mCI = './mu_mCI'
filename_mu_pCI = './mu_pCI'
filename_p_0CI = './p_0CI'
filename_tauCI = './tauCI'
filename_sigma1CI = './sigma1CI'
filename_sigma2CI = './sigma2CI'

filename_curve1_upper = './curve1_upper'
filename_curve2_upper = './curve2_upper'
filename_curve1_lower = './curve1_lower'
filename_curve2_lower = './curve2_lower'

filename_iter = './iterations'

ifAppend=FALSE;
for(iter in 1:NITER){
  if(iter == 1){
    ifAppend=FALSE
  }else{
    ifAppend=TRUE
  }
  set.seed(seeds[iter])
  source('mainDDE.R')
  
  write.table(matrix(round(meanmu_m, digits = 5),nr=1), file = paste(filename_mu_m,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meanmu_p, digits = 5),nr=1), file = paste(filename_mu_p,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meanp_0, digits = 5),nr=1), file = paste(filename_p_0,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meantau, digits = 5),nr=1), file = paste(filename_tau,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meansigma1, digits = 5),nr=1), file = paste(filename_sigma1,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meansigma2, digits = 5),nr=1), file = paste(filename_sigma2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mu_mCredible, digits = 5),nr=1), file = paste(filename_mu_mCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mu_pCredible, digits = 5),nr=1), file = paste(filename_mu_pCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(p_0Credible, digits = 5),nr=1), file = paste(filename_p_0CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(tauCredible, digits = 5),nr=1), file = paste(filename_tauCI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);  
  write.table(matrix(round(sigma1Credible, digits = 5),nr=1), file = paste(filename_sigma1CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(sigma2Credible, digits = 5),nr=1), file = paste(filename_sigma2CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mean_curve1, digits = 5),nr=1), file = paste(filename_curve1,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(mean_curve2, digits = 5),nr=1), file = paste(filename_curve2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(iter, digits = 5),nr=1), file = paste(filename_iter,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(lambda, digits = 5),nr=1), file = paste(filename_lambda,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
  write.table(matrix(round(upper_curve1, digits = 5),nr=1), file = paste(filename_curve1_upper,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(upper_curve2, digits = 5),nr=1), file = paste(filename_curve2_upper,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(lower_curve1, digits = 5),nr=1), file = paste(filename_curve1_lower,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(lower_curve2, digits = 5),nr=1), file = paste(filename_curve2_lower,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
}


