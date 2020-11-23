#rm(list=ls())
#library("fda")
#library("MCMCpack")
#library(numDeriv)
#library(fda)
#library(deSolve)
#library(NLRoot)
#library(smcUtils)
#library(Hmisc)

#times <- seq(0, 100, 0.2)
times <- seq(0, 500, 5)
sigma1 <- 1
sigma2 <- 10
sigma1_simulator <- sigma1
sigma2_simulator <- sigma2
DDE_power <- 8
source("splineDDE.R")
source('simulatorDDE.R')

set.seed(364)
M <- nknots
alambda <- 1
blambda <- 1
sigmac2 <- 1
#CESSthresholds <- 0.8
#NP <- 300

mu_m_particles = abs(rnorm(NP, 0.04, 0.02))
mu_p_particles = abs(rnorm(NP, 0.04, 0.02))
p_0_particles = rnorm(NP, 90, 20)
tau_particles = rnorm(NP, 30, 10)

#sigma1 <- sqrt(1/rgamma(NP, 1,1))
#sigma2 <- sqrt(1/rgamma(NP, 1,1))

sigma1 <- runif(NP, 0.3, 3)
sigma2 <- runif(NP, 0.2, 8)

lambda <- (rgamma(NP, 1, 1))

c1_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
c2_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2

#sum((output[,2]-basismat%*%c1_LSE)^2)
double_quadpts <- table(quadpts)[which(table(quadpts) == 2)]

##optimize ###
objectiveFunction <- function(par, lambda){
  mu_m = par[1]
  mu_p = par[2]
  p_0 = par[3]
  tau = par[4]
  c1 <- par[5:((length(par))/2+2)]
  c2 <- par[-(1:((length(par))/2+2))]
  rt1 <- sum((y1 - basismat%*%c1)^2/0.3^2+(y2 - basismat%*%c2)^2/3^2)
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  index <- min(which(quadpts >= tau))
  if(quadpts[index+1] > quadpts[index]){
    x2_delay <- c(rep(x2[1], index), x2[2:(nrow(D0quadbasismat)-index+1)])
  }else{
    x2_delay <- c(rep(x2[1], index+1), x2[2:(nrow(D0quadbasismat)-index)])
  }
  g1 <- 1/(1+(x2_delay/p_0)^10) - mu_m*x1 
  g2 <- x1 - mu_p*x2
  rt2 <- lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  rt <- rt1 + rt2
  return(rt)
}

#par <- c(0.04, 0.04, 90, 25.0, c1_LSE, c2_LSE)
#optimizer <- nlminb(start = par, objectiveFunction, control=list(trace=FALSE, abs.tol = 0.1, iter.max = 100), lambda = 1)
#c1_LSE <- optimizer$par[3:((length(par))/2+1)]
#c2_LSE <- optimizer$par[-(1:((length(par))/2+1))]

c1 <- matrix(NA, nr = NP, nc = ncol(basismat))
c2 <- matrix(NA, nr = NP, nc = ncol(basismat))

for(i in 1:NP){
  c1[i,] <- rnorm(length(c1_LSE), c1_LSE, 0.1)
  c2[i,] <- rnorm(length(c2_LSE), c2_LSE, 0.1)
}

#resampleThreshold <- 0.5

source('functionDDE.R')
source('asmcDDE.R')
