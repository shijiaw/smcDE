rm(list=ls())
library(deSolve)
#library("fda")
library(Hmisc)
#library("MCMCpack")
library(numDeriv)
library(fda)
library(NLRoot)
library(smcUtils)
library(MASS)

NITER <- 40
knotsNpoints <- 6
#lambda <- 10
NP <- 1000
lambda <- rep(1, NP)
CESSthresholds <- 0.95
resampleThreshold <- 0.5
seeds <- 1:NITER
filename_kappa1 = './kappa1'
filename_kappa2 = './kappa2'
filename_sigma1 = './sigma1'
filename_sigma2 = './sigma2'
filename_curve1 = './curve1'
filename_curve2 = './curve2'
filename_lambda = './lambda'

filename_kappa1CI = './kappa1CI'
filename_kappa2CI = './kappa2CI'
filename_sigma1CI = './sigma1CI'
filename_sigma2CI = './sigma2CI'

filename_iter = './iterations'

ifAppend=FALSE;
for(iter in 1:NITER){
  if(iter == 1){
    ifAppend=FALSE
  }else{
    ifAppend=TRUE
  }
  set.seed(seeds[iter])
  source('mainBasis.R')
  
  write.table(matrix(round(meankappa1, digits = 5),nr=1), file = paste(filename_kappa1,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meankappa2, digits = 5),nr=1), file = paste(filename_kappa2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meansigma1, digits = 5),nr=1), file = paste(filename_sigma1,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(meansigma2, digits = 5),nr=1), file = paste(filename_sigma2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(kappa1Credible, digits = 5),nr=1), file = paste(filename_kappa1CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(kappa2Credible, digits = 5),nr=1), file = paste(filename_kappa2CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(sigma1Credible, digits = 5),nr=1), file = paste(filename_sigma1CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(sigma2Credible, digits = 5),nr=1), file = paste(filename_sigma2CI,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(MeanSpline1, digits = 5),nr=1), file = paste(filename_curve1,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(MeanSpline2, digits = 5),nr=1), file = paste(filename_curve2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(iter, digits = 5),nr=1), file = paste(filename_iter,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  write.table(matrix(round(lambda, digits = 5),nr=1), file = paste(filename_lambda,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
}


