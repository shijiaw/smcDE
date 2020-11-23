## Simulator ###
#set.seed(452)
#library(deSolve)

parameters <- c(r = 0.8, P = 2, tau = 3)
parameters

initial <- c(3500)
initial
#state
DE_model <- function(t, y, parameters) {
  ## assign parameter values to parameter variables
  ## the function 'unname()' just removes the name of the parameter - this is unnecessary, it just cleans things up a bit
  r <- unname(parameters['r'])
  P <- unname(parameters['P'])
  tau <- unname(parameters['tau'])
  
  ## define time lags ###
  if (t < tau)
    lag <- initial
  else
    lag <- lagvalue(t - tau)
  
  ## compute the rates of change
  #dX1dt <- 1/(1+(lag/p_0)^2) - mu_m*X1   
  #dX2dt <- X1 - mu_p*X2
     
  #dy1 <- 1/(1+(lag/p_0)^DDE_power) - mu_m*y[1]
  dy <- r*y*(1-lag/(1000*P))
  
  return(list(dy, dy = dy))
  #return(list(c(dy)))

}

#times <- seq(0, 200, 1)

## Simulate the system!
output <- dede(y=initial, times=times, func=DE_model, parms=parameters, atol = 1e-10)


## Simulate observations y ##
y <- exp(as.vector(log(output[,2]) + rnorm(length(times), 0, sigma1)))

# gname = c("DDEobs.eps",sep="")  
# postscript(gname,width=10,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
# par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
# 
# plot(output[,1:2], lty = 1, lwd = 2, type = 'l', xlab = expression(X[1](t)), ylab = '',ylim = c(0, 10000))
# points(output[,1], y1, col = 2, cex = 0.2)
# 
# dev.off()


# theta <- c(0.03, 0.03, 100, 25)
# ####ODE w.r.t parameters ###
# DDEf <- function(theta){
#   parameters <- c(mu_m = theta[1], mu_p = theta[2], p_0 = theta[3], tau = theta[4])
#   #parameters <- c(kappa1 = theta[1], kappa2 = theta[2])
#   output <- dede(y=initial, times=times, func=DE_model, parms=parameters)
#   return(list(X1 = output[,2], X2 = output[,3]))
# }
# 
# ###Likelihood function w.r.t theta ###
# l <- function(theta, y1, y2, sigma1, sigma2){
#   #theta <- c(theta1, theta2)
#   output <- DDEf(theta)
#   X1 <- as.vector(output$X1)
#   X2 <- as.vector(output$X2)
#   sum <- -sum((y1 - X1)^2/(2*sigma1^2) + (y2 - X2)^2/(2*sigma2^2))
#   return(sum)
# }
# 
# #l(theta, y1, y2, 1, 3)
# l_sigma1 <- function(theta, y1, y2, sigma1, sigma2){
#   #theta <- c(theta1, theta2)
#   output <- DDEf(theta)
#   X1 <- as.vector(output$X1)
#   X2 <- as.vector(output$X2)
#   sum <- -sum((y1 - X1)^2/(2*sigma1^2)) - length(y1)*log(sigma1)
#   return(sum)
# }

#l(c(0.030, 0.030, 100, 25), y1, y2, 1, 3)
#l_sigma1(c(0.030, 0.030, 100, 25), y1, y2, 1, 3)
