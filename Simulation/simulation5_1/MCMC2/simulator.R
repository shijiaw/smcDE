## Simulator ###
set.seed(452)
library(deSolve)

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
  dX1dt <- 72/(36+X2) - abs(kappa1)
  dX2dt <- kappa2*X1 - 1
  
  ## return as a list object, with the first element of the list being the derivatives. The order of derivatives must be the same as the order in the initial condition vector!
  return(list(c(dX1dt, dX2dt)))
}

times <- seq(0, 60, 0.5)

## Simulate the system!
output <- ode(y=initial, times=times, func=DE_model, parms=parameters)


## Simulate observations y ##
y1 <- as.vector(output[,2] + rnorm(121, 0, 1))
y2 <- as.vector(output[,3] + rnorm(121, 0, 3))

####ODE w.r.t parameters ###
ODEf <- function(theta, initialvalue){
  parameters <- c(kappa1 = theta[1], kappa2 = theta[2])
  initial <- c(X1=initialvalue[1], X2=initialvalue[2])
  output <- ode(y=initial, times=times, func=DE_model, parms=parameters, method = "euler")
  return(list(X1 = output[,2], X2 = output[,3]))
}

###Likelihood function w.r.t theta ###
l <- function(theta, y1, y2, sigma1, sigma2){
  #theta <- c(theta1, theta2)
  output <- ODEf(theta)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  sum <- -sum((y1 - X1)^2/sigma1^2 + (y2 - X2)^2/sigma2^2)
  return(sum)
}

tempLfun <- function(theta, initial, sigma, y1, y2){
  temperature <- 1
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  output <- ODEf(theta, initial)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  sum <- -sum((y1 - X1)^2/sigma1^2 + (y2 - X2)^2/sigma2^2)*temperature-length(y1)*log(sigma1^2)*temperature-length(y2)*log(sigma2^2)*temperature
  return(sum)
}

###Metropolis Hastings algorithm ###
MHTheta <- function(oldTheta, Initial, y1, y2, sigma){
  temperature = 1
  MHsigma1 <- 0.012
  newTheta <- c(rnorm(1, oldTheta[1], MHsigma1), rnorm(1, oldTheta[2], MHsigma1))
  #newInitial <- c(rnorm(1, oldInitial[1], MHsigma2), rnorm(1, oldInitial[2], MHsigma2))
  logacc_prob <- tempLfun(newTheta, Initial, sigma, y1, y2) - tempLfun(oldTheta, Initial, sigma, y1, y2)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- newTheta
  }else{
    rt <- oldTheta
  }
  return(rt)
}

MHInitial <- function(Theta, oldInitial, y1, y2, sigma){
  MHsigma2 <- 0.05
  newInitial <- c(rnorm(1, oldInitial[1], MHsigma2), rnorm(1, oldInitial[2], MHsigma2))
  logacc_prob <- tempLfun(Theta, newInitial, sigma, y1, y2) - tempLfun(Theta, oldInitial, sigma, y1, y2)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- newInitial
  }else{
    rt <- oldInitial
  }
  return(rt)
}

##Gibbs move for sigma ###
GBKernel <- function(theta, initial, y1, y2){
  output <- ODEf(theta, initial)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  temp1 <- rgamma(1, shape = 1+ length(y1)/2, rate = 1 + sum((y1 - X1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2, rate = 1 + sum((y2 - X2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}




