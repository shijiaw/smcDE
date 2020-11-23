library(MASS)

Loglikelihood <- function(theta, sigma, y1, y2, c1, c2){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]

  x1 <- basismat%*%c1
  x2 <- basismat%*%c2
  
  sum <- -sum((y1 - x1)^2/sigma1^2 + (y2 - x2)^2/sigma2^2) - length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2)
  return(sum)
}




tempLfun <- function(theta, sigma, y1, y2, c1, c2){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  
  x1 <- basismat%*%c1
  x2 <- basismat%*%c2
  
  sum <- -sum((y1 - x1)^2/sigma1^2 + (y2 - x2)^2/sigma2^2) -length(y1)*log(sigma1^2) - length(y2)*log(sigma2^2) 
  return(sum)
}

tempPosterior <- function(kappa, sigma, y1, y2, c1, c2, lambda){
  sigma1 <- as.numeric(sigma[1])
  sigma2 <- as.numeric(sigma[2])
  
  xx1 <- basismat%*%c1
  xx2 <- basismat%*%c2
  
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - abs(kappa[1]) 
  g2 <- kappa[2]*x1 - 1
  
  rt <- -sum((y1 - xx1)^2/sigma1^2 + (y2 - xx2)^2/sigma2^2) -length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2) - lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}


###Metropolis Hastings algorithm ###
MH_kappa <- function(lambda, kappaOld, sigma, c1, c2){
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1Old <- 72/(36+x2) - abs(kappaOld[1]) 
  g2Old <- kappaOld[2]*x1 - 1
  
  # u <- runif(1)
  # if(u < 0.25){
  #   MHsigma <- 0.02
  # }
  # if(u >= 0.25 && u < 0.5){
  #   MHsigma <- 0.06
  # }
  # if(u >= 0.5 && u < 0.75){
  #   MHsigma <- 0.2
  # }
  # if(u >= 0.75){
  #   MHsigma <- 0.6
  # }
  MHsigma <- 0.2
  kappanew <- c(rnorm(1, kappaOld[1], MHsigma), rnorm(1, kappaOld[2], MHsigma))
  
  g1new <- 72/(36+x2) - abs(kappanew[1]) 
  g2new <- kappanew[2]*x1 - 1
  
  logacc_prob <- (-lambda*sum(quadwts*((dx1-g1new)^2+(dx2-g2new)^2)) + lambda*sum(quadwts*((dx1-g1Old)^2+(dx2-g2Old)^2)))

  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(kappanew)
  }else{
    rt <- c(kappaOld)
  }
  return(rt)
}


c2_ll <- function(c2, y2, sigma2){
  rt <- -sum((y2 - basismat%*%c2)^2/sigma2^2)
  return(rt)
}

c2_logprior <- function(lambda, kappa, c1, c2, y2){
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - abs(kappa[1]) 
  g2 <- kappa[2]*x1 - 1
  rt <- -lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}

MH_c2 <- function(c2_priorOld, c2_llOld, lambda, sigma2, kappa, c1, c2, y2){
  MHsigma <- 1
  newc2 <- c(rnorm(length(c2), c2, sigmac2))
  c2_priornew <- c2_logprior(lambda, kappa, c1, c2, y2)
  c2_llnew <- c2_ll(newc2, y2, sigma2)
  logacc_prob <- (c2_priornew - c2_priorOld) + (c2_llnew - c2_llOld)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(newc2)
  }else{
    rt <- c(c2)
  }
  return(list(c2_prior = c2_priornew, c2_ll = c2_llnew, c2 =rt)) 
}

##Gibbs move for sigma ###
GB_sigma <- function(y1, y2, c1, c2){
  x1 <- basismat%*%c1
  x2 <- basismat%*%c2
  temp1 <- rgamma(1, shape = 1+ length(y1)/2, rate = 1 + sum((y1 - x1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2, rate = 1 + sum((y2 - x2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}

##Gibbs move for lambda ###
GB_lambda <- function(alambda, blambda, kappa, M, c1, c2){
  kappa1 <- kappa[1]
  kappa2 <- kappa[2]
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - abs(kappa1) 
  g2 <- kappa2*x1 - 1
  b_lambda <- 1/(1/blambda+sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))/2)
  return(rgamma(1, shape = alambda+M/2, scale = b_lambda))
}


##Gibbs move for c1 ###
GB_c1 <- function(lambda, sigma1, kappa, c1, c2, y1){
  VV <- diag(quadwts)
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - abs(kappa[1]) 
  g2 <- kappa[2]*x2 - 1
  SigInv <- (lambda/2*(kappa[2]^2*t(D0quadbasismat)%*%VV%*%D0quadbasismat+t(D1quadbasismat)%*%VV%*%D1quadbasismat)+ t(basismat)%*%basismat/(2*sigma1^2))
  Sig <- solve(SigInv)
  mu <- t((lambda/2*(t(dx2+1)%*%VV%*%D0quadbasismat*kappa[2]+t(g1)%*%VV%*%D1quadbasismat)+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
  return(mvrnorm(n = 1, mu, Sig))
}


