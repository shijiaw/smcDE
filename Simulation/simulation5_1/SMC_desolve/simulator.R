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
  dX1dt <- 72/(36+X2) - kappa1
  dX2dt <- kappa2*X1 - 1
  
  ## return as a list object, with the first element of the list being the derivatives. The order of derivatives must be the same as the order in the initial condition vector!
  return(list(c(dX1dt, dX2dt)))
}

times <- seq(0, 60, 0.5)

## Simulate the system!
output <- ode(y=initial, times=times, func=DE_model, parms=parameters)
#plot(output, lty = 1, lwd = 2)

## Simulate observations y ##
y1 <- as.vector(output[,2] + rnorm(121, 0, 1))
y2 <- as.vector(output[,3] + rnorm(121, 0, 3))


