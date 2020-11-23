rm(list=ls())
library("fda")
#library("MCMCpack")
library(numDeriv)
library(fda)
library(deSolve)
times <- seq(0.5, 60, 0.5)
source("spline_basis.R")

times <- c(0, times)
###
parameters <- c(kappa1 = 2, kappa2 = 1)
parameters

initial <- c(X1=7, X2=-10)
initial

DE_model <- function(t, state, parameters) {
  ## assign parameter values to parameter variables
  ## the function 'unname()' just removes the name of the parameter - this is unnecessary, it just cleans things up a bit
  kappa1 <- unname(parameters['kappa1'])
  kappa2 <- unname(parameters['kappa2'])
  
  ## assign state variables to names
  X1 <- unname(state['X1'])
  X2 <- unname(state['X2'])
  
  ## compute the rates of change
  dX1dt <- 72/(36+X2) - kappa1
  dX2dt <- kappa2*X1 - 1
  
  ## return as a list object, with the first element of the list being the derivatives. The order of derivatives must be the same as the order in the initial condition vector!
  return(list(c(dX1dt, dX2dt)))
}


## Simulate the system!
output <- ode(y=initial, times=times, func=DE_model, parms=parameters)
#plot(output, lty = 1, lwd = 2)

## Simulate observations y ##
y1 <- as.vector(output[-1,2] + rnorm(120, 0, 1))
y2 <- as.vector(output[-1,3] + rnorm(120, 0, 2))

c1_initial <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
c2_initial <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2
fitted <- basismat%*%c1_initial
plot(y1, type = 'l')
points(fitted, col = 2) 

##optimize ###
objectiveFunction <- function(par, lambda){
  kappa1 <- par[1]
  kappa2 <- par[2]
  c1 <- par[3:((length(par))/2+1)]
  c2 <- par[-(1:((length(par))/2+1))]
  rt1 <- sum((y1 - basismat%*%c1)^2/1^2+(y2 - basismat%*%c2)^2/2^2)
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - kappa1 
  g2 <- kappa2*x1 - 1
  rt2 <- lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  rt <- rt1 + rt2
  return(rt)
}

par <- c(5, 7, c1_initial, c2_initial)
optimizer <- nlminb(start = par, objectiveFunction, control=list(trace=FALSE, abs.tol = 0.1, iter.max = 100), lambda = 1000)
#gamma2_temp <- rt$par
kappa1_estimate <- optimizer$par[1]
kappa2_estimate <- optimizer$par[2]
c1_estimate <- optimizer$par[3:((length(par))/2+1)]
c2_estimate <- optimizer$par[-(1:((length(par))/2+1))]

fitted1 <- basismat%*%c1_estimate
fitted2 <- basismat%*%c2_estimate
plot(output[-1,2], type = 'l')
points(fitted1, col = 2)

plot(output[-1,3], type = 'l')
points(fitted2, col = 2)


# lambda0 = 1
# for(iter in 1:10){
#   lambda = lambda*5
#   par <- c(gamma_estimate, c1_estimate, c2_estimate)
#   optimizer <- nlminb(start = par, objectiveFunction, control=list(trace=FALSE, abs.tol = 0.1, iter.max = 100), lambda = lambda)
#   #gamma2_temp <- rt$par
#   gamma_estimate <- optimizer$par[1]
#   print(gamma_estimate)
#   c1_estimate <- optimizer$par[2:((length(par)+1)/2)]
#   c2_estimate <- optimizer$par[-(1:((length(par)+1)/2))]
# }
# ###posterior surface plot



