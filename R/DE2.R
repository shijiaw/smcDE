DE2 <- function(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac){
  y = data[[1]]
  times = times[[1]]
  
  # set random seed for smc 
  set.seed(seed)
  
  tobs = times
  knots = knots
  nknots   = length(knots)
  norder   = 4
  nbasis   = length(knots) + norder - 2
  bsbasis = create.bspline.basis(range(knots), nbasis, norder, knots)
  
  # basis values at sampling points
  basismat = eval.basis(tobs, bsbasis)
  
  # square of basis matrix (Phi). It is used frequent for optimization,
  # so we calculate in advance
  basismat2  = t(basismat)%*%basismat
  
  # values of the first derivative of basis functions at sampling points
  Dbasismat  = eval.basis(tobs, bsbasis, 1)
  
  # values of the second derivative of basis functions at sampling points
  D2basismat = eval.basis(tobs, bsbasis, 2)
  
  # number of quadrature points (including two knots) between two knots
  nquad = 5;
  
  # set up quadrature points and weights
  #[sinbasis, quadpts, quadwts]
  quadre = quadset(nquad, bsbasis)
  quadpts=quadre$quadvals[,1]
  quadwts=quadre$quadvals[,2]
  
  
  # values of the second derivative of basis functions at quadrature points
  D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
  D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
  
  
  psi <- basismat  # bsplineS(norder=4,nderiv=0)
  psi1 <- Dbasismat # bsplineS(norder=4,nderiv=1)
  
  M <- nknots
  
  # initialize particles 
  # Note that the particles can also be generated from other distributions (e.g. priors)
  r_particles = runif(NP, 0, 2)
  P_particles = runif(NP, 0, 5)
  tau_particles = runif(NP, 0, 10)
  
  sigma <- runif(NP, 0.3, 5)
  
  lambda <- rgamma(NP, 1, 10)
  
  c_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%log(y)
  
  c <- matrix(NA, nr = NP, nc = ncol(basismat))
  
  for(i in 1:NP){
    c[i,] <- rnorm(length(c_LSE), c_LSE, 0.2)
  }
  
  # adaptive annealing SMC, use B-spline to solve ODE
  ESSVec <- list()
  parametersStore <- list()
  sigmaStore <- list()
  tempSeq <- list()
  tempDifSeq <- list()
  cStore <- list()
  lambdaStore <- list()
  
  
  logZ = 0
  alpha <- 0
  alphaDif <- 0
  
  particles <- cbind(c, r_particles, P_particles, tau_particles, sigma, lambda)
  
  # weights function 
  W <- rep(1/NP, NP)
  logW <- log(W)
  iter <- 0
  
  while(alpha < 1){
    iter = iter + 1
    cat("iteration:", iter, "\n")
    
    # Compute weights and conduct resampling 
    c <- particles[,1:ncol(basismat)] 
    sigma <- particles[,ncol(basismat)+4] 
    parameters <- particles[,(ncol(basismat)+1):(ncol(basismat)+3)]
    
    logRatio <- rep(NA, NP)
   
    for(i in 1:NP){
      logRatio[i] <- tempPosterior(as.vector(parameters[i,]), sigma[i], y, as.vector(c[i,]), lambda[i], basismat, knots, bsbasis, nquad)
    }
    
    alphaDif <- bisection(0, 1, W, logRatio, CESSthresholds)
    alpha <- alpha + alphaDif
    cat("annealing parameter:", alpha, "\n")
    if(alpha > 1){
      alpha <- 1
      alphaDif <- 1- tempSeq[[iter-1]]
    }
    tempSeq[[iter]] <- alpha
    tempDifSeq[[iter]] <- alphaDif
    
    logw <- alphaDif*logRatio + logW
    logmax = max(logw)
    W = exp(logw-logmax)/sum(exp(logw-logmax)) 
    logW = log(W)
    logZ = logZ + log(sum(exp(logw - logmax))) + logmax 
    RESS <- ESS(logW)
    ESSVec[[iter]] <- RESS
    
    if(RESS < resampleThreshold){
      ancestor <- systematic.resample(W, num.samples = length(W), engine="R")
      particesOld = particles[ancestor,]
      logW = rep(-log(NP), NP)
      W = rep(1/NP, NP)
    }else{
      particesOld = particles
    }
    
    ptm <- proc.time()
    
    c <- particles[, 1:ncol(basismat)]
    parameters <- particles[, (ncol(basismat)+1):(ncol(basismat)+3)]
    sigma <- particles[, (ncol(basismat)+4)]
    
    for(i in 1:NP){
      # Evaluate Basis functions 
      Basis_eval <- newBasis(knots, parameters[i,3], bsbasis, nquad)
      quadpts_eval <- Basis_eval$quadpts
      quadwts_eval <- Basis_eval$quadwts
      D0quadbasismat_eval <- Basis_eval$D0quadbasismat
      D1quadbasismat_eval <- Basis_eval$D1quadbasismat
      D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
      
      # Evaluate the curves 
      ww <- basismat%*%c[i,]
      dw <- D1quadbasismat_eval%*%c[i,]
      w <- D0quadbasismat_eval%*%c[i,]
      tau_index = which(quadpts_eval == parameters[i,3])[1]
      w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
      
      # sample lambda using Gibbs
      lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, w, dw, w_delay, quadwts_eval, alpha)
      
      # sample sigma
      sigma[i] <- GB_sigma(y, ww, alpha)
      
      # sample c
      c_llOld <- c_ll(ww, y, sigma[i])
      c_priorOld <- c_logprior(lambda[i], parameters[i,], dw, w_delay, quadwts_eval)
      MHstep_c <- MH_c(c_priorOld, c_llOld, lambda[i], sigma[i], parameters[i,], c[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y, alpha)
      c[i,] <- MHstep_c$c
      
      ww <- basismat%*%c[i,]
      dw <- D1quadbasismat_eval%*%c[i,]
      w <- D0quadbasismat_eval%*%c[i,]
      tau_index = which(quadpts_eval == parameters[i,3])[1]
      w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
      
      # sample parameters
      parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i], c[i,], dw, w, w_delay, quadwts_eval, alpha, knots, bsbasis, nquad)
    }
    
    proc.time() - ptm
    
    parametersStore[[iter]] <- parameters
    sigmaStore[[iter]] <- sigma
    cStore[[iter]] <- c
    lambdaStore[[iter]] <- lambda
    
    particles <- cbind(c, parameters, sigma, lambda)
  }
  
  # particles at final iteration
  c <- particles[, 1:ncol(basismat)]
  parameters <- particles[, (ncol(basismat)+1):(ncol(basismat)+3)]
  sigma <- particles[, (ncol(basismat)+4)]
  lambda <- particles[, (ncol(basismat)+5)]
  
  parameterList <- list(r = parameters[,1], P = parameters[,2], tau = parameters[,3])
  
  
  return(list(c = c, parameters = parameterList, lambda = lambda, sigma = sigma, W = W))
}

# function for evaluating basis functions with new tau sampled
newBasis <- function(knots, tau, bsbasis, nquad){
  new_eval = sort(c(knots, tau))
  quadre = quadset(nquad, bsbasis, new_eval)
  quadpts = quadre$quadvals[,1]
  
  tau_index = which(quadpts == tau)[1]
  
  D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
  D0quadbasismat_delay = eval.basis((quadpts-tau)[-(1:tau_index)], bsbasis, 0)
  D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
  return(list(quadpts=quadpts, quadwts=quadre$quadvals[,2], D0quadbasismat = D0quadbasismat, D1quadbasismat = D1quadbasismat, D0quadbasismat_delay = D0quadbasismat_delay))
}

# Compute tempering posterior
tempPosterior <- function(parameters, sigma, y, c, lambda, basismat, knots, bsbasis, nquad){
  y <- log(y)
  
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  
  sigma <- as.numeric(sigma)
  
  Basis_eval <- newBasis(knots, tau, bsbasis, nquad)
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


# Metropolis Hastings sampling for parameters (r, p, tau)
MH_parameter <- function(lambda, parametersOld, sigma, c, dw, w, w_delay, quadwts, temperature, knots, bsbasis, nquad){
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
  # taunew <- rnorm(1, tauOld, tausigma)
  taunew <- rtruncnorm(1, a=0, b=10, mean = tauOld, sd = tausigma)
  
  ## re-define basis eval ##
  Basis_eval <- newBasis(knots, taunew, bsbasis, nquad)
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

# log likelihood for c
c_ll <- function(ww, y, sigma){
  y <- log(y)
  #w <- basismat%*%c
  rt <- -sum((y - ww)^2/sigma^2)
  return(rt)
}

# prior distribution for c
c_logprior <- function(lambda, parameters, dw, w_delay, quadwts){
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  
  g <- r*(1-exp(w_delay)/(1000*P))
  
  rt <- -lambda*sum(quadwts*((dw-g)^2))
  return(rt)
}

# Metropolis Hasting to sample basis coefficients c
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

# Gibbs step for measurement error sigma 
GB_sigma <- function(y, ww, temperature){
  y <- log(y)
  temp <- rgamma(1, shape = 1+ length(y)/2, rate = 1 + temperature*sum((y - ww)^2)/2)
  return(sqrt(1/temp))
}

# Gibbs move for lambda 
GB_lambda <- function(alambda, blambda, parameters, M, w, dw,  w_delay, quadwts, temperature){
  r <- parameters[1]
  P <- parameters[2]
  tau <- parameters[3]
  g <- r*(1-exp(w_delay)/(1000*P))
  b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dw-g)^2))/2)
  return(rgamma(1, shape = alambda+temperature*M/2, scale = b_lambda))
}


# effective sample size
ESS <- function(logW){
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(NP)
  return(exp(logRESS))
}


bisection <- function(a, b, W, logRatio, alphaCESS) {
  n = 1000
  tol = 1e-7
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if(func(W, logRatio, a, alphaCESS)*func(W, logRatio, b, alphaCESS) > 0){
    stop('signs of f(a) and f(b) differ')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((func(W, logRatio, c, alphaCESS) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(func(W, logRatio, b, alphaCESS)) == sign(func(W, logRatio, a, alphaCESS)), 
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
  rt <- (sum(w*W))^2/(sum(W*w^2))- alphaCESS
  return(rt)
}


