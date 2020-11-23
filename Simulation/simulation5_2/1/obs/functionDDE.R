library(MASS)

### Basis Evaluation ###
newBasis <- function(knots, tau, bsbasis){
  new_eval = sort(c(knots, tau))
  quadre = quadset(nquad, bsbasis, new_eval)
  quadpts = quadre$quadvals[,1]
  
  tau_index = which(quadpts == tau)[1]
  
  D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
  D0quadbasismat_delay = eval.basis((quadpts-tau)[-(1:tau_index)], bsbasis, 0)
  D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
  return(list(quadpts=quadpts, quadwts=quadre$quadvals[,2], D0quadbasismat = D0quadbasismat, D1quadbasismat = D1quadbasismat, D0quadbasismat_delay = D0quadbasismat_delay))
}


tempPosterior <- function(parameters, sigma, y, c, lambda){
  y <- log(y)
  
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  
  sigma <- as.numeric(sigma)
  
  Basis_eval <- newBasis(knots, tau, bsbasis)
  quadpts_eval <- Basis_eval$quadpts
  quadwts_eval <- Basis_eval$quadwts
  D0quadbasismat_eval <- Basis_eval$D0quadbasismat
  D1quadbasismat_eval <- Basis_eval$D1quadbasismat
  D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay

  ww <- basismat%*%c
  
  dw <- D1quadbasismat_eval%*%c
  w <- D0quadbasismat_eval%*%c
  
  tau_index = which(quadpts_eval == tau)[1]
  w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c)
  
  g <- r*(1-exp(w_delay)/(1000*P))
  
  rt <- -sum((y - ww)^2/sigma^2) -length(y)*log(sigma^2) - lambda*sum(quadwts_eval*((dw-g)^2))
  return(rt)
}


###Metropolis Hastings algorithm ###
MH_parameter <- function(lambda, parametersOld, sigma, c, dw, w, w_delay, quadwts, temperature){
  r_Old <- parametersOld[1]
  P_Old <- parametersOld[2]
  tauOld <- parametersOld[3]
  Q_Old <- r_Old/P_Old
  
  sigmar_Inv = temperature*lambda*sum(quadwts)
  sigmar = 1/sqrt(sigmar_Inv)
  mu_r = temperature*lambda*(sum(quadwts*(dw+Q_Old*exp(w_delay)/1000)))/sigmar_Inv
  rnew <- rtruncnorm(1, a=0, b=Inf, mean = mu_r, sd = sigmar)

  sigmaQ_Inv <- temperature*lambda*sum(quadwts*(exp(w_delay)/(1000))^2)
  sigmaQ <- 1/sqrt(sigmaQ_Inv)
  mu_Q <- -temperature*lambda*(sum(quadwts*(dw-rnew)*(exp(w_delay)/(1000))))/sigmaQ_Inv
  Qnew <- rtruncnorm(1, a=0, b=Inf, mean = mu_Q, sd = sigmaQ)
  P_new <- rnew/Qnew

  g1Old <- rnew*(1-exp(w_delay)/(1000*P_new))
  
  u2 <- runif(1)
  if(u2 < 0.5){
    tausigma <- 0.2
  }
  if(u2 >= 0.5){
    tausigma <- 2
  }
  #taunew <- rnorm(1, tauOld, tausigma)
  taunew <- rtruncnorm(1, a=0, b=10, mean = tauOld, sd = tausigma)
  
  ## re-define basis eval ##
  Basis_eval <- newBasis(knots, taunew, bsbasis)
  quadpts_eval <- Basis_eval$quadpts
  quadwts_eval <- Basis_eval$quadwts
  D0quadbasismat_eval <- Basis_eval$D0quadbasismat
  D1quadbasismat_eval <- Basis_eval$D1quadbasismat
  D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
  
  dw_new <- D1quadbasismat_eval%*%c
  w_new <- D0quadbasismat_eval%*%c
  tau_index = which(quadpts_eval == taunew)[1]
  w_delay <- c(rep(w_new[1], tau_index), D0quadbasismat_delay%*%c)
  
  g1new <- rnew*(1-exp(w_delay)/(1000*P_new))
  
  logacc_prob <- temperature*(-lambda*sum(quadwts_eval*((dw_new-g1new)^2)) + lambda*sum(quadwts*((dw-g1Old)^2)))
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(rnew, P_new, taunew)
  }else{
    rt <- c(rnew, P_new, tauOld)
  }  
  
  return(rt)
}


c_ll <- function(ww, y, sigma){
  y <- log(y)
  #w <- basismat%*%c
  rt <- -sum((y - ww)^2/sigma^2)
  return(rt)
}

c_logprior <- function(lambda, parameters, dw, w_delay, quadwts){
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  
  # w <- D0quadbasismat%*%c
  # dw <- D1quadbasismat%*%c
  # 
  # index <- min(which(quadpts >= tau))
  # if(quadpts[index+1] > quadpts[index]){
  #   w_delay <- c(rep(w[1], index), w[2:(nrow(D0quadbasismat)-index+1)])
  # }else{
  #   w_delay <- c(rep(w[1], index+1), w[2:(nrow(D0quadbasismat)-index)])
  # }
  g <- r*(1-exp(w_delay)/(1000*P))
  
  rt <- -lambda*sum(quadwts*((dw-g)^2))
  return(rt)
}

MH_c <- function(c_priorOld, c_llOld, lambda, sigma, parameters, c, quadwts, tau_index, basismat, D0quadbasismat, D1quadbasismat, D0quadbasismat_delay, y, temperature){
  #MHsigma <- 1
  newc <- c(rnorm(length(c), c, sigmac))
  
  ww <- basismat%*%newc
  dw <- D1quadbasismat%*%newc
  w <- D0quadbasismat%*%newc
  w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%newc)
  
  c_priornew <- c_logprior(lambda, parameters, dw, w_delay, quadwts)
  c_llnew <- c_ll(ww, y, sigma)
  logacc_prob <- temperature*(c_priornew - c_priorOld) + temperature*(c_llnew - c_llOld)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    #print("accept")
    rt <- c(newc)
  }else{
    #print("reject")
    rt <- c(c)
  }
  return(list(c_prior = c_priornew, c_ll = c_llnew, c =rt)) 
}

##Gibbs move for sigma ###
GB_sigma <- function(y, ww, temperature){
  y <- log(y)
  #w <- basismat%*%c
  temp <- rgamma(1, shape = 1+ length(y)/2, rate = 1 + temperature*sum((y - ww)^2)/2)
  return(sqrt(1/temp))
}

##Gibbs move for lambda ###
GB_lambda <- function(alambda, blambda, parameters, M, w, dw,  w_delay, quadwts, temperature){
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  
  # w <- D0quadbasismat%*%c
  # dw <- D1quadbasismat%*%c
  # 
  # index <- min(which(quadpts >= tau))
  # if(quadpts[index+1] > quadpts[index]){
  #   w_delay <- c(rep(w[1], index), w[2:(nrow(D0quadbasismat)-index+1)])
  # }else{
  #   w_delay <- c(rep(w[1], index+1), w[2:(nrow(D0quadbasismat)-index)])
  # }
  g <- r*(1-exp(w_delay)/(1000*P))
  
  #rt <- -lambda*sum(quadwts*((dw-g)^2))
  #g1 <- 72/(36+x2) - abs(kappa1) 
  #g2 <- kappa2*x1 - 1
  b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dw-g)^2))/2)
  return(rgamma(1, shape = alambda+temperature*M/2, scale = b_lambda))
}


# ##Gibbs move for c1 ###
# GB_c1 <- function(lambda, sigma1, parameters, c1, c2, y1, temperature){
#   mu_m <- parameters[1]
#   mu_p <- parameters[2]
#   p_0 <- parameters[3]
#   tau <- parameters[4]
#   
#   VV <- diag(quadwts)
#   dx1 <- D1quadbasismat%*%c1
#   dx2 <- D1quadbasismat%*%c2
#   x1 <- D0quadbasismat%*%c1
#   x2 <- D0quadbasismat%*%c2
#   
#   index <- min(which(quadpts >= tau))
#   if(quadpts[index+1] > quadpts[index]){
#     x2_delay <- c(rep(x2[1], index), x2[2:(nrow(D0quadbasismat)-index+1)])
#   }else{
#     x2_delay <- c(rep(x2[1], index+1), x2[2:(nrow(D0quadbasismat)-index)])
#   }
#   #g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
#   #g2 <- x1 - mu_p*x2
# 
#   SigInv <- temperature*(lambda/2*(t(D0quadbasismat)%*%VV%*%D0quadbasismat+t(D1quadbasismat+mu_m*D0quadbasismat)%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+ t(basismat)%*%basismat/(2*sigma1^2))
#   Sig <- solve(SigInv)
#   mu <- t(temperature*(lambda/2*(t(dx2+mu_p*x2)%*%VV%*%D0quadbasismat+t(1/(1+(x2_delay/p_0)^DDE_power))%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
#   return(mvrnorm(n = 1, mu, Sig))
# }

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

