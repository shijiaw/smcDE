
Loglikelihood <- function(theta, sigma, y1, y2, c1, c2){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]

  x1 <- basismat%*%c1
  x2 <- basismat%*%c2
  
  sum <- -sum((y1 - x1)^2/sigma1^2 + (y2 - x2)^2/sigma2^2) - length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2)
  return(sum)
}




tempLfun <- function(theta, sigma, y1, y2, c1, c2, temperature){
  sigma1 <- sigma[1]
  sigma2 <- sigma[2]
  
  x1 <- basismat%*%c1
  x2 <- basismat%*%c2
  
  sum <- -sum((y1 - x1)^2/sigma1^2 + (y2 - x2)^2/sigma2^2)*temperature-length(y1)*log(sigma1^2)*temperature-length(y2)*log(sigma2^2)*temperature
  return(sum)
}

tempPosterior <- function(kappa, sigma, y1, y2, xx1, xx2, dx1, dx2, x1, x2, lambda){
  sigma1 <- as.numeric(sigma[1])
  sigma2 <- as.numeric(sigma[2])
  
  # xx1 <- basismat%*%c1
  # xx2 <- basismat%*%c2
  # 
  # dx1 <- D1quadbasismat%*%c1
  # dx2 <- D1quadbasismat%*%c2
  # x1 <- D0quadbasismat%*%c1
  # x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - kappa[1] 
  g2 <- kappa[2]*x1 - 1
  
  rt <- -sum((y1 - xx1)^2/sigma1^2 + (y2 - xx2)^2/sigma2^2) -length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2) - lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}


###Metropolis Hastings algorithm ###
#Gibbs_kappa <- function(lambda, sigma, c1, c2, temperature){
Gibbs_kappa <- function(lambda, sigma, dx1, dx2, x1, x2, temperature){
  # dx1 <- D1quadbasismat%*%c1
  # dx2 <- D1quadbasismat%*%c2
  # x1 <- D0quadbasismat%*%c1
  # x2 <- D0quadbasismat%*%c2

  sigma1_Inv = lambda*sum(quadwts*1)
  sigma1 = 1/sqrt(sigma1_Inv)
  mu_kappa1 = lambda*(sum(quadwts*(72/(36+x2)-dx1)))/sigma1_Inv
  kappa1new <- rnorm(1, mean = mu_kappa1, sd = sigma1)
  
  sigma2_Inv = lambda*sum(quadwts*x1^2)
  sigma2 = 1/sqrt(sigma2_Inv)
  mu_kappa2 = lambda*(sum(quadwts*(x1*(dx2+1))))/sigma2_Inv
  kappa2new <- rnorm(1, mean = mu_kappa2, sd = sigma2)
  
  return(c(kappa1new, kappa2new))
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
  g1 <- 72/(36+x2) - kappa[1]
  g2 <- kappa[2]*x1 - 1
  rt <- -lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}

MH_c2 <- function(c2_priorOld, c2_llOld, lambda, sigma2, kappa, c1, c2, y2, temperature){
  MHsigma <- 2
  newc2 <- c(rnorm(length(c2), c2, sigmac2))
  c2_priornew <- c2_logprior(lambda, kappa, c1, c2, y2)
  c2_llnew <- c2_ll(newc2, y2, sigma2)
  logacc_prob <- temperature*(c2_priornew - c2_priorOld) + temperature*(c2_llnew - c2_llOld)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(newc2)
    #print('accept')
  }else{
    rt <- c(c2)
    #print('reject')
  }
  return(list(c2_prior = c2_priornew, c2_ll = c2_llnew, c2 =rt)) 
}

##Gibbs move for sigma ###
GB_sigma <- function(y1, y2, xx1, xx2, temperature){
  # x1 <- basismat%*%c1
  # x2 <- basismat%*%c2
  temp1 <- rgamma(1, shape = 1+ length(y1)/2*temperature, rate = 1 + temperature*sum((y1 - xx1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2*temperature, rate = 1 + temperature*sum((y2 - xx2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}

##Gibbs move for lambda ###
GB_lambda <- function(alambda, blambda, kappa, M, dx1, dx2, x1, x2, temperature){
  kappa1 <- kappa[1]
  kappa2 <- kappa[2]
  g1 <- 72/(36+x2) - kappa1
  g2 <- kappa2*x1 - 1
  b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))/2)
  return(rgamma(1, shape = alambda+temperature*M*2/2, scale = b_lambda))
}


##Gibbs move for c1 ###
GB_c1 <- function(lambda, sigma1, kappa, dx2, x2, y1, temperature){
  g1 <- 72/(36+x2) - kappa[1]
  #g2 <- kappa[2]*x1 - 1
  SigInv <- temperature*(lambda/2*(kappa[2]^2*D0VVD0+D1VVD1)+ bsmat2/(2*sigma1^2))
  Sig <- solve(SigInv)
  mu <- t(temperature*(lambda/2*(t(dx2+1)%*%VV%*%D0quadbasismat*kappa[2]+t(g1)%*%VV%*%D1quadbasismat)+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
  return(mvrnorm(n = 1, mu, Sig))
}

ESS <- function(logW){
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(NP)
  return(exp(logRESS))
  #return(1/sum(w^2)/length(w))
}

# CESS <- function(w, W){
#   return(length(w)*(sum(w*W))^2/(sum(W*w^2)))
# }



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


