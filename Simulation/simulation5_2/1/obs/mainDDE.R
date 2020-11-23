rm(list=ls())
library("fda")
library("truncnorm")
library(numDeriv)
library(fda)
library(deSolve)
library(NLRoot)
library(Hmisc)


times <- seq(0, 100, 0.25)
knots = seq(min(times), max(times),by=2)
sigma1 <- 0.4
sigma1_simulator <- sigma1

set.seed(123)
source("splineDDE.R")
source('simulatorDDE.R')

set.seed(364)
M <- nknots
alambda <- 1
blambda <- 1
sigmac <- 0.3
CESSthresholds <- 0.3
NP <- 500


r_particles = abs(rnorm(NP, 1.2, 1))
P_particles = abs(rnorm(NP, 2, 1))
tau_particles = rnorm(NP, 2, 0.6)


sigma <- runif(NP, 0.3, 5)

lambda <- (rgamma(NP, 1, 1))

c_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%log(y)



double_quadpts <- table(quadpts)[which(table(quadpts) == 2)]



c <- matrix(NA, nr = NP, nc = ncol(basismat))

for(i in 1:NP){
  c[i,] <- rnorm(length(c_LSE), c_LSE, 0.2)
}

resampleThreshold <- 0.5

source('functionDDE.R')
source('resample.R')
source('asmcDDE.R')
