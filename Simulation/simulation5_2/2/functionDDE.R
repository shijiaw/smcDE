library(MASS)


### Basis Evaluation ###
newBasis <- function(knots, tau, bsbasis){
  new_eval = sort(c(knots, tau))
  quadre = quadset(nquad, bsbasis, new_eval)
  quadpts=quadre$quadvals[,1]
  
  tau_index = which(quadpts == tau)[1]
  
  D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
  D0quadbasismat_delay = eval.basis((quadpts-tau)[-(1:tau_index)], bsbasis, 0)
  D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
  return(list(quadpts=quadpts, quadwts=quadre$quadvals[,2], D0quadbasismat = D0quadbasismat, D1quadbasismat = D1quadbasismat, D0quadbasismat_delay = D0quadbasismat_delay))
}


tempPosterior <- function(parameters, sigma, y1, y2, c1, c2, lambda){
  mu_m <- parameters[1]
  mu_p <- parameters[2]
  p_0 <- parameters[3]
  tau <- parameters[4]
  
  sigma1 <- as.numeric(sigma[1])
  sigma2 <- as.numeric(sigma[2])
  
  Basis_eval <- newBasis(knots, tau, bsbasis)
  quadpts_eval <- Basis_eval$quadpts
  quadwts_eval <- Basis_eval$quadwts
  D0quadbasismat_eval <- Basis_eval$D0quadbasismat
  D1quadbasismat_eval <- Basis_eval$D1quadbasismat
  D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
  
  xx1 <- basismat%*%c1
  xx2 <- basismat%*%c2
  dx1 <- D1quadbasismat_eval%*%c1
  dx2 <- D1quadbasismat_eval%*%c2
  x1 <- D0quadbasismat_eval%*%c1
  x2 <- D0quadbasismat_eval%*%c2
  tau_index = which(quadpts_eval == tau)[1]
  x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2)
  
  g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
  g2 <- x1 - mu_p*x2
  
  rt <- -sum((y1 - xx1)^2/sigma1^2 + (y2 - xx2)^2/sigma2^2) -length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2) - lambda*sum(quadwts_eval*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}


###Metropolis Hastings algorithm ###
MH_parameter <- function(lambda, parametersOld, sigma, c1, c2, dx1, dx2, x1, x2, x2_delay, quadwts, temperature){
  mu_mOld <- parametersOld[1]
  mu_pOld <- parametersOld[2]
  p_0Old <- parametersOld[3]
  tauOld <- parametersOld[4]
  
  temp_delay <- (x2_delay/p_0Old)^DDE_power
  
  sigmamu_mInv = temperature*lambda*sum(quadwts*x1^2)
  sigmamu_m = 1/sqrt(sigmamu_mInv)
  mu_mu_m = -temperature*lambda*(sum(quadwts*x1*(dx1-1/(1+temp_delay))))/sigmamu_mInv
  mu_mnew <- rtruncnorm(1, a=0, b=100, mean = mu_mu_m, sd = sigmamu_m)
  
  sigmamu_pInv = temperature*lambda*sum(quadwts*x2^2)
  sigmamu_p = 1/sqrt(sigmamu_pInv)
  mu_mu_p = temperature*lambda*(sum(quadwts*x2*(x1-dx2)))/sigmamu_pInv
  mu_pnew <- rtruncnorm(1, a=0, b=100, mean = mu_mu_p, sd = sigmamu_p)
  
  u1 <- runif(1)
  if(u1 < 0.5){
    psigma <- 2
  }
  if(u1 >= 0.5){
    psigma <- 10
  }
  
  p_0new <- rtruncnorm(1, a=0, b=Inf, mean = p_0Old, sd = psigma)
    
  g1Old <- 1/(1+ (x2_delay/p_0Old)^DDE_power) - mu_mnew*x1 
  g1new <- 1/(1+ (x2_delay/p_0new)^DDE_power) - mu_mnew*x1 
  #g2new <- x1 - mu_pnew*x2
  
  logacc_prob <- dnorm(p_0new, 90, 20, log = TRUE) - dnorm(p_0Old, 90, 20, log = TRUE) +temperature*(-lambda*sum(quadwts*((dx1-g1new)^2)) + lambda*sum(quadwts*((dx1-g1Old)^2)))
  #logacc_prob <- tempLfun(newTheta, newInitial, sigma, y1, y2, temperature) - tempLfun(oldTheta, oldInitial, sigma, y1, y2, temperature)
  if(runif(1, 0, 1) > exp(logacc_prob)){
    p_0new <- p_0Old
  }
  
  g1Old <- 1/(1+(x2_delay/p_0new)^DDE_power) - mu_mnew*x1 
  
  u2 <- runif(1)
  if(u2 < 0.5){
    tausigma <- 1
  }
  if(u2 >= 0.5){
    tausigma <- 5
  }
  taunew <- rtruncnorm(1, a=0, b=50, mean = tauOld, sd = tausigma)
  
  ## re-define basis eval ##
  Basis_eval <- newBasis(knots, taunew, bsbasis)
  quadpts_eval <- Basis_eval$quadpts
  quadwts_eval <- Basis_eval$quadwts
  D0quadbasismat_eval <- Basis_eval$D0quadbasismat
  D1quadbasismat_eval <- Basis_eval$D1quadbasismat
  D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
  
  dx1_new <- D1quadbasismat_eval%*%c1
  dx2_new <- D1quadbasismat_eval%*%c2
  x1_new <- D0quadbasismat_eval%*%c1
  x2_new <- D0quadbasismat_eval%*%c2
  tau_index = which(quadpts_eval == taunew)[1]
  x2_delay <- c(rep(x2_new[1], tau_index), D0quadbasismat_delay%*%c2)
  
  g1new <- 1/(1+(x2_delay/p_0new)^DDE_power) - mu_mnew*x1_new 
  #g2new <- x1 - mu_pnew*x2
  
  
  logacc_prob <- temperature*(-lambda*sum(quadwts_eval*((dx1_new-g1new)^2)) + lambda*sum(quadwts*((dx1-g1Old)^2)))
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(mu_mnew, mu_pnew, p_0new, taunew)
  }else{
    rt <- c(mu_mnew, mu_pnew, p_0new, tauOld)
  }  
  
  return(rt)
}


c2_ll <- function(xx2, y2, sigma2){
  rt <- -sum((y2 - xx2)^2/sigma2^2)
  return(rt)
}


c2_logprior <- function(lambda, parameters, dx1, dx2, x1, x2, x2_delay, quadwts, y2){
  mu_m <- parameters[1]
  mu_p <- parameters[2]
  p_0 <- parameters[3]
  tau <- parameters[4]
  
  g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
  g2 <- x1 - mu_p*x2
  
  rt <- -lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}

MH_c2 <- function(c2_priorOld, c2_llOld, lambda, sigma2, parameters, dx1, x1, c2, quadwts, tau_index, basismat, D0quadbasismat, D1quadbasismat, D0quadbasismat_delay, y2, temperature){
  MHsigma <- 1
  newc2 <- c(rnorm(length(c2), c2, sigmac2))
  
  xx2 <- basismat%*%newc2
  dx2 <- D1quadbasismat%*%newc2
  x2 <- D0quadbasismat%*%newc2
  x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%newc2)
  
  c2_priornew <- c2_logprior(lambda, parameters, dx1, dx2, x1, x2, x2_delay, quadwts, y2)
  c2_llnew <- c2_ll(xx2, y2, sigma2)
  logacc_prob <- temperature*(c2_priornew - c2_priorOld) + temperature*(c2_llnew - c2_llOld)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(newc2)
  }else{
    rt <- c(c2)
  }
  return(list(c2_prior = c2_priornew, c2_ll = c2_llnew, c2 =rt)) 
}

##Gibbs move for sigma ###
GB_sigma <- function(y1, y2, x1, x2, temperature){
  temp1 <- rgamma(1, shape = 1+ length(y1)/2, rate = 1 + temperature*sum((y1 - x1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2, rate = 1 + temperature*sum((y2 - x2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}

##Gibbs move for lambda ###
GB_lambda <- function(alambda, blambda, parameters, M, dx1, dx2, x1, x2,  x2_delay, quadwts, temperature){
  mu_m <- parameters[1]
  mu_p <- parameters[2]
  p_0 <- parameters[3]
  tau <- parameters[4]
  
  g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
  g2 <- x1 - mu_p*x2
  
  b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))/2)
  return(rgamma(1, shape = alambda+temperature*M*2/2, scale = b_lambda))
}


##Gibbs move for c1 ###
GB_c1 <- function(lambda, sigma1, parameters, dx1, dx2, x1, x2,  x2_delay, quadwts, D0quadbasismat, D1quadbasismat, y1, temperature){
  mu_m <- parameters[1]
  mu_p <- parameters[2]
  p_0 <- parameters[3]
  tau <- parameters[4]
  
  VV <- diag(quadwts)

  SigInv <- temperature*(lambda/2*(t(D0quadbasismat)%*%VV%*%D0quadbasismat+t(D1quadbasismat+mu_m*D0quadbasismat)%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+ t(basismat)%*%basismat/(2*sigma1^2))
  Sig <- solve(SigInv)
  mu <- t(temperature*(lambda/2*(t(dx2+mu_p*x2)%*%VV%*%D0quadbasismat+t(1/(1+(x2_delay/p_0)^DDE_power))%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
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


