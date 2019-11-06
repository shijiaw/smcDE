DE1 <- function(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac){
  y1 = data[[1]]
  y2 = data[[2]]
  times = times[[1]]
  sigmac2 = sigmac
  # set random seed for smc 
  set.seed(seed)
  
  tobs = times
  knots = knots
  nknots   = length(knots)
  norder   = 4
  nbasis   = length(knots) + norder - 2
  bsbasis = create.bspline.basis(range(knots),nbasis,norder,knots)
  
  # basis values at sampling points
  basismat   = eval.basis(tobs, bsbasis)
  
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
  
  VV <- diag(quadwts)
  D0VVD0 <- t(D0quadbasismat)%*%VV%*%D0quadbasismat
  D1VVD1 <- t(D1quadbasismat)%*%VV%*%D1quadbasismat
  bsmat2 <- t(basismat)%*%basismat
  
  M <- nknots
  
  theta1 <- rnorm(NP, 4, 3)
  theta2 <- rnorm(NP, 4, 3)
  
  sigma1 <- sqrt(1/rgamma(NP, 1,1))
  sigma2 <- sqrt(1/rgamma(NP, 1,1))
  
  
  c1_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
  c2_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2
  
  double_quadpts <- table(quadpts)[which(table(quadpts) == 2)]
  
  c1 <- matrix(NA, nr = NP, nc = ncol(basismat))
  c2 <- matrix(NA, nr = NP, nc = ncol(basismat))
  
  for(i in 1:NP){
    c1[i,] <- rnorm(length(c1_LSE), c1_LSE, 1)
    c2[i,] <- rnorm(length(c2_LSE), c2_LSE, 3)
  }
  
  lambda <- rgamma(NP, 1, 1)
  
  # adaptive annealing SMC, use B-spline to solve ODE 
  ESSVec <- list()
  
  kappaStore <- list()
  sigmaStore <- list()
  tempSeq <- list()
  tempDifSeq <- list()
  c1Store <- list()
  c2Store <- list()

  dx1 <- matrix(NA, nr = nrow(D1quadbasismat), nc = NP)
  x1 <- matrix(NA, nr = nrow(D1quadbasismat), nc = NP)
  xx1 <- matrix(NA, nr = nrow(basismat), nc = NP)
  
  # x2, dx2, xx2 
  dx2 <- matrix(NA, nr = nrow(D1quadbasismat), nc = NP)
  x2 <- matrix(NA, nr = nrow(D1quadbasismat), nc = NP)
  xx2 <- matrix(NA, nr = nrow(basismat), nc = NP)
  
  logZ = 0
  alpha <- 0
  alphaDif <- 0
  particles <- cbind(c1, c2, theta1, theta2, sigma1, sigma2)
  
  W <- rep(1/NP, NP)
  logW <- log(W)
  iter <- 0

  while(alpha < 1){
    iter = iter + 1
    cat("iteration:", iter, "\n")
    # Compute weights and conduct resampling 
    c1 <- particles[,1:ncol(basismat)] 
    c2 <- particles[,(ncol(basismat)+1):(2*ncol(basismat))] 
    
    for(i in 1:NP){
      dx1[,i] <- D1quadbasismat%*%c1[i,]
      x1[,i] <- D0quadbasismat%*%c1[i,]
      xx1[,i] <- basismat%*%c1[i,]
      
      dx2[,i] <- D1quadbasismat%*%c2[i,]
      x2[,i] <- D0quadbasismat%*%c2[i,]
      xx2[,i] <- basismat%*%c2[i,]
      
    }
    
    sigma1 <- particles[,2*ncol(basismat)+3] 
    sigma2 <- particles[,2*ncol(basismat)+4]
    sigma <- cbind(sigma1, sigma2)
    kappa <- particles[,(2*ncol(basismat)+1):(2*ncol(basismat)+2)]
    logRatio <- rep(NA, NP)
    for(i in 1:NP){
      logRatio[i] <- tempPosterior(as.vector(kappa[i,]), as.vector(sigma[i,]), y1, y2, xx1[,i], xx2[,i], dx1[,i], dx2[,i], x1[,i], x2[,i], lambda[i], quadwts)
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
    
    kappa <- particles[, (2*ncol(basismat)+1):(2*ncol(basismat)+2)]
    sigma1 <- particles[, (2*ncol(basismat)+3)]
    sigma2 <- particles[, (2*ncol(basismat)+4)]
    sigma <- cbind(sigma1, sigma2)
    for(i in 1:NP){
      # sample lambda using Gibbs
      lambda[i] <- GB_lambda(alambda, blambda, kappa[i,], M, dx1[,i], dx2[,i], x1[,i], x2[,i], alpha, quadwts)
      
      # sample sigma
      sigma[i,] <- GB_sigma(y1, y2, xx1[,i], xx2[,i], alpha)
      
      # sample kappa
      kappa[i,] <- Gibbs_kappa(lambda[i], sigma[i,], dx1[,i], dx2[,i], x1[,i], x2[,i], alpha, quadwts)
      
      # sample c1
      c1[i,] <- GB_c1(lambda[i], sigma1[i], kappa[i,], dx2[,i], x2[,i], y1, alpha, basismat,D0quadbasismat, D1quadbasismat, D0VVD0, D1VVD1, VV, bsmat2)
      
      # sample c2
      c2_llOld <- c2_ll(c2[i,], y2, sigma2[i], basismat)
      c2_priorOld <- c2_logprior(lambda[i], kappa[i,], c1[i,], c2[i,], y2, D0quadbasismat, D1quadbasismat, quadwts)
      MHstep_c2 <- MH_c2(c2_priorOld, c2_llOld, lambda[i], sigma2[i], kappa[i,], c1[i,], c2[i,], y2, alpha, D0quadbasismat, D1quadbasismat, quadwts, basismat)
      c2[i,] <- MHstep_c2$c2
      
    }
    proc.time() - ptm
    
    particles <- cbind(c1, c2, kappa, sigma1, sigma2)
    
  }
  parameters <- list(kappa1 = kappa[,1], kappa2 = kappa[,2])
  
  return(list(c1 = c1, c2 = c2, parameters = parameters, lambda = lambda, sigma = sigma, W = W))
}


tempPosterior <- function(kappa, sigma, y1, y2, xx1, xx2, dx1, dx2, x1, x2, lambda, quadwts){
  sigma1 <- as.numeric(sigma[1])
  sigma2 <- as.numeric(sigma[2])

  g1 <- 72/(36+x2) - kappa[1] 
  g2 <- kappa[2]*x1 - 1
  
  rt <- -sum((y1 - xx1)^2/sigma1^2 + (y2 - xx2)^2/sigma2^2) -length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2) - lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}


# Metropolis Hastings algorithm 
Gibbs_kappa <- function(lambda, sigma, dx1, dx2, x1, x2, temperature, quadwts){
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


c2_ll <- function(c2, y2, sigma2, basismat){
  rt <- -sum((y2 - basismat%*%c2)^2/sigma2^2)
  return(rt)
}

c2_logprior <- function(lambda, kappa, c1, c2, y2, D0quadbasismat, D1quadbasismat, quadwts){
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - kappa[1]
  g2 <- kappa[2]*x1 - 1
  rt <- -lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt)
}

MH_c2 <- function(c2_priorOld, c2_llOld, lambda, sigma2, kappa, c1, c2, y2, temperature, D0quadbasismat, D1quadbasismat, quadwts, basismat){
  MHsigma <- 1
  newc2 <- c(rnorm(length(c2), c2, MHsigma))
  c2_priornew <- c2_logprior(lambda, kappa, c1, c2, y2, D0quadbasismat, D1quadbasismat, quadwts)
  c2_llnew <- c2_ll(newc2, y2, sigma2, basismat)
  logacc_prob <- temperature*(c2_priornew - c2_priorOld) + temperature*(c2_llnew - c2_llOld)
  if(runif(1, 0, 1) < exp(logacc_prob)){
    rt <- c(newc2)
  }else{
    rt <- c(c2)
  }
  return(list(c2_prior = c2_priornew, c2_ll = c2_llnew, c2 =rt)) 
}

# Gibbs move for sigma 
GB_sigma <- function(y1, y2, xx1, xx2, temperature){
  temp1 <- rgamma(1, shape = 1+ length(y1)/2*temperature, rate = 1 + temperature*sum((y1 - xx1)^2)/2)
  temp2 <- rgamma(1, shape = 1+ length(y2)/2*temperature, rate = 1 + temperature*sum((y2 - xx2)^2)/2)
  return(c(sqrt(1/temp1), sqrt(1/temp2)))
}

# Gibbs move for lambda 
GB_lambda <- function(alambda, blambda, kappa, M, dx1, dx2, x1, x2, temperature, quadwts){
  kappa1 <- kappa[1]
  kappa2 <- kappa[2]
  g1 <- 72/(36+x2) - kappa1
  g2 <- kappa2*x1 - 1
  b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))/2)
  return(rgamma(1, shape = alambda+temperature*M/2, scale = b_lambda))
}


# Gibbs move for c1 
GB_c1 <- function(lambda, sigma1, kappa, dx2, x2, y1, temperature,basismat, D0quadbasismat, D1quadbasismat, D0VVD0, D1VVD1, VV, bsmat2){
  g1 <- 72/(36+x2) - kappa[1]
  SigInv <- temperature*(lambda/2*(kappa[2]^2*D0VVD0+D1VVD1)+ bsmat2/(2*sigma1^2))
  Sig <- solve(SigInv)
  mu <- t(temperature*(lambda/2*(t(dx2+1)%*%VV%*%D0quadbasismat*kappa[2]+t(g1)%*%VV%*%D1quadbasismat)+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
  return(mvrnorm(n = 1, mu, Sig))
}

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


