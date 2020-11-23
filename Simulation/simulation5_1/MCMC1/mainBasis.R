rm(list=ls())
library("fda")
library(numDeriv)
library(fda)
library(deSolve)
library(NLRoot)

times <- seq(0, 60, 0.5)
source("spline_basis.R")
source('simulator.R')

set.seed(364)
M <- nknots
alambda <- 1
blambda <- 1
sigmac2 <- 0.1
Niter <- 400000
kappa <- c(4,3)


sigma <- c(2, 6)


lambda <- (1/rgamma(1, 1, 10))

c1_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
c2_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2

c1 <- rnorm(length(c1_LSE), rep(0, length(c1_LSE)), 3)
c2 <- rnorm(length(c2_LSE), rep(0, length(c2_LSE)), 5)


kappaStore <- matrix(NA, nr = Niter, nc = 2)
sigmaStore <- matrix(NA, nr = Niter, nc = 2)
c1Store <- matrix(NA, nr = Niter, nc = 18)
c2Store <- matrix(NA, nr = Niter, nc = 18)
lambdaStore <- rep(NA, Niter)

kappaStore[1,] <- kappa
sigmaStore[1,] <- sigma
c1Store[1,] <- c1
c2Store[1,] <- c2
lambdaStore[1] <- lambda

source('functionBasisExpansion.R')
source('MHBasis.R')

mean(kappaStore[-(1:200000), 1])
mean(kappaStore[-(1:200000), 2])
quantile(kappaStore[-(1:200000), 1], c(0.025, 0.975))
quantile(kappaStore[-(1:200000), 2], c(0.025, 0.975))

mean(c1Store[-(1:200000), 1])
mean(c2Store[-(1:200000), 1])
quantile(c1Store[-(1:200000), 1], c(0.025, 0.975))
quantile(c2Store[-(1:200000), 1], c(0.025, 0.975))

