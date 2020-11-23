#rm(list=ls())

times <- seq(0, 60, 0.5)
source("spline_basis.R")
source('simulator.R')

M <- nknots
alambda <- 1
blambda <- 1
sigmac2 <- 0.1

theta1 <- rnorm(NP, 4, 3)
theta2 <- rnorm(NP, 4, 3)

sigma1 <- sqrt(1/rgamma(NP, 1,1))
sigma2 <- sqrt(1/rgamma(NP, 1,1))


c1_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
c2_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2


c1 <- matrix(NA, nr = NP, nc = ncol(basismat))
c2 <- matrix(NA, nr = NP, nc = ncol(basismat))


for(i in 1:NP){
  c1[i,] <- rnorm(length(c1_LSE), c1_LSE, 1)
  c2[i,] <- rnorm(length(c2_LSE), c2_LSE, 3)
}


source('functionBasisExpansion.R')
source('asmcBasis.R')
