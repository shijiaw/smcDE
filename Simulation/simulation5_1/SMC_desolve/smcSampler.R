### SMC_Sampler ###
library(deSolve)
library(Hmisc)
times <- seq(0, 60, 0.5)
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

####ODE w.r.t parameters ###
ODEf <- function(theta, initial){
  Initial <- c(X1=initial[1], X2=initial[2])
  parameters <- c(kappa1 = theta[1], kappa2 = theta[2])
  output <- ode(y=Initial, times=times, func=DE_model, parms=parameters, method = "euler")
  return(list(X1 = output[,2], X2 = output[,3]))
}


tempLfun <- function(theta, initial, sigma, y1, y2, temperature){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  output <- ODEf(theta, initial)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  sum <- -sum((y1 - X1)^2/sigma1^2 + (y2 - X2)^2/sigma2^2)*temperature-length(y1)*log(sigma1^2)*temperature-length(y2)*log(sigma2^2)*temperature
  return(sum)
}

Lfun <- function(theta, initial, sigma, y1, y2){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  output <- ODEf(theta, initial)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  sum <- -sum((y1 - X1)^2/sigma1^2 + (y2 - X2)^2/sigma2^2) -length(y1)*log(sigma1^2) - length(y2)*log(sigma2^2)
  return(sum)
}

bisection <- function(f, a, b, W, logRatio, alphaCESS) {
  n = 1000
  tol = 1e-7
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if(f(W, logRatio, a, alphaCESS)*f(W, logRatio, b, alphaCESS) > 0){
    stop('signs of f(a) and f(b) differ')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(W, logRatio, c, alphaCESS) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(W, logRatio, b, alphaCESS)) == sign(f(W, logRatio, a, alphaCESS)), 
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print('Too many iterations')
}

func <- function(W, logRatio, a, alphaCESS) {
  logw <- a*logRatio
  logmax = max(logw)
  w <- exp(logw - logmax)
  #rt <- length(w)*(sum(w*W))^2/(sum(W*w^2))- alphaCESS
  rt <- (sum(w*W))^2/(sum(W*w^2))- alphaCESS
  return(rt)
}

###Metropolis Hastings algorithm ###
MHKernel <- function(oldTheta, oldInitial, y1, y2, sigma, temperature){
  #sigma1 <- sigma[1]
  #sigma1 <- sigma[2]
  MHsigma1 <- 0.04
  MHsigma2 <- 0.5
  newTheta <- c(rnorm(1, oldTheta[1], MHsigma1), rnorm(1, oldTheta[2], MHsigma1))
  newInitial <- c(rnorm(1, oldInitial[1], MHsigma2), rnorm(1, oldInitial[2], MHsigma2))
  logacc_prob <- tempLfun(newTheta, newInitial, sigma, y1, y2, temperature) - tempLfun(oldTheta, oldInitial, sigma, y1, y2, temperature)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(newTheta, newInitial, sigma)
  }else{
    rt <- c(oldTheta, oldInitial, sigma)
  }
  return(rt)
}

##Gibbs move for sigma ###
GBKernel <- function(theta, initial, y1, y2, temperature){
  output <- ODEf(theta, initial)
  X1 <- as.vector(output$X1)
  X2 <- as.vector(output$X2)
  temp1 <- rgamma(1, shape = 1+ length(y1)/2*temperature, rate = 1 + temperature*sum((y1 - X1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2*temperature, rate = 1 + temperature*sum((y2 - X2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}

set.seed(364)
NP <- 500
CESSthresholds <- 0.999
theta1 <- runif(NP, -4, 4)
theta2 <- runif(NP, 0, 4)

sigma1 <- runif(NP, 0.3, 6)
sigma2 <- runif(NP, 0.3, 6)

initial1 <- runif(NP, -15, 15)
initial2 <- runif(NP, -15, 15)


ESS <- function(w){
  return(1/sum(w^2)/length(w))
}

ESSVec <- rep(NA, T)

thetaStore <- list()
sigmaStore <- list()
initialStore <- list()

logZ = 0

resampleThreshold <- 0.5
particles <- cbind(theta1, theta2, initial1, initial2, sigma1, sigma2)
logW = -log(NP)
W = rep(1/NP, NP)

iter = 0
alpha <- 0
alphaDif <- 0
tempSeq <- list()
tempDifSeq <- list()

#for(t in tempScheme){
while(alpha < 1){
  print(alpha)
  iter = iter + 1
  RESS <- ESS(W)
  ESSVec[iter] <- RESS
  if(RESS < resampleThreshold){
    ancestor = sample.int(NP, prob = W, replace = TRUE)
    particesOld = particles[ancestor,]
    logW = -log(NP)
    W = rep(1/NP, NP)
  }else{
    particesOld = particles
  }
  
  logRatio <- rep(NA, NP)
  logRatio <- apply(as.matrix(particles), 1, function(x)(Lfun(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))
  #for(i in 1:NP){
  #  logRatio[i] <- apply(as.matrix(particles), 1, function(x)(Lfun(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))
  #}
  alphaDif <- bisection(func, 0, 1, W, logRatio, CESSthresholds)
  
  alpha <- alpha + alphaDif
  
  if(alpha > 1){
    alpha <- 1
    alphaDif <- 1- tempSeq[[iter-1]]
  }
  tempSeq[[iter]] <- alpha
  tempDifSeq[[iter]] <- alphaDif
  
  ##Reweighting
  if(alpha < 1){
    #logw <- apply(as.matrix(particles), 1, function(x)(tempLfun(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2, tempDif[iter+1])))+logW
    logw <- alphaDif*logRatio + logW
    logmax = max(logw)
    logZ = logZ + log(sum(exp(logw - logmax))) + logmax 
    
    W = exp(logw-logmax)/sum(exp(logw-logmax)) 
    RESS <- ESS(W)
    logW = log(W)
    
  }

  ##MH move
  particles <- t(apply(particesOld, 1, function(x)(MHKernel(as.vector(x)[1:2], as.vector(x)[3:4], y1, y2, as.vector(x)[5:6], alpha))))
  particles[,5:6] <- t(apply(particesOld, 1, function(x)(GBKernel(as.vector(x)[1:2], as.vector(x)[3:4], y1, y2, alpha))))
  
  thetaStore[[iter]] <- particles[,1:2]
  sigmaStore[[iter]] <- particles[,5:6]
  initialStore[[iter]] <- particles[,3:4]
  apply(thetaStore[[iter]], 2, mean)
  apply(sigmaStore[[iter]], 2, mean)
  apply(initialStore[[iter]], 2, mean)
  #particleStore[[t]] <- particles
  
}


meantheta1 <- sum(abs(thetaStore[[iter]][,1])*W)
meantheta2 <- sum(thetaStore[[iter]][,2]*W)
meansigma1 <- sum(sigmaStore[[iter]][,1]*W)
meansigma2 <- sum(sigmaStore[[iter]][,2]*W)
meaninitial1 <- sum(initialStore[[iter]][,1]*W)
meaninitial2 <- sum(initialStore[[iter]][,2]*W)
print(meantheta1)
print(meantheta2)
print(meansigma1)
print(meansigma2)
print(meaninitial1)
print(meaninitial2)


theta1Credible <- wtd.quantile(abs(thetaStore[[iter]][,1]), na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
theta2Credible <- wtd.quantile(thetaStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
sigma1Credible <- wtd.quantile(sigmaStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
sigma2Credible <- wtd.quantile(sigmaStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
initial1Credible <- wtd.quantile(initialStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
initial2Credible <- wtd.quantile(initialStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))




