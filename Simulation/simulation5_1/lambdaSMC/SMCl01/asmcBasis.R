### adaptive annealing SMC, use B-spline to solve ODE###
ESSVec <- list()

kappaStore <- list()
sigmaStore <- list()
tempSeq <- list()
tempDifSeq <- list()
c1Store <- list()
c2Store <- list()


# x1, dx1, xx1 
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
  ### Compute weights and conduct resampling ###
  #logRatio <- apply(as.matrix(particles), 1, function(x)(Loglikelihood(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))
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
    logRatio[i] <- tempPosterior(as.vector(kappa[i,]), as.vector(sigma[i,]), y1, y2, xx1[,i], xx2[,i], dx1[,i], dx2[,i], x1[,i], x2[,i], lambda[i])
  }
  alphaDif <- bisection(func, 0, 1, W, logRatio, CESSthresholds)
  
  alpha <- alpha + alphaDif
  
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
    
    ##sample sigma
    sigma[i,] <- GB_sigma(y1, y2, xx1[,i], xx2[,i], alpha)
    
    ##sample kappa
    kappa[i,] <- Gibbs_kappa(lambda[i], sigma[i,], dx1[,i], dx2[,i], x1[,i], x2[,i], alpha)
    #kappa[i,] <- MH_kappa(lambda, kappa[i,], sigma[i,], c1[i,], c2[i,], alpha)
    
    ##sample c1
    c1[i,] <- GB_c1(lambda[i], sigma1[i], kappa[i,], dx2[,i], x2[,i], y1, alpha)
    
    ##sample c2
    c2_llOld <- c2_ll(c2[i,], y2, sigma2[i])
    c2_priorOld <- c2_logprior(lambda[i], kappa[i,], c1[i,], c2[i,], y2)
    MHstep_c2 <- MH_c2(c2_priorOld, c2_llOld, lambda[i], sigma2[i], kappa[i,], c1[i,], c2[i,], y2, alpha)
    c2[i,] <- MHstep_c2$c2
    
  }
  proc.time() - ptm
  
  kappaStore[[iter]] <- kappa
  sigmaStore[[iter]] <- sigma
  c1Store[[iter]] <- c1
  c2Store[[iter]] <- c2
  #lambdaStore[[iter]] <- lambda
  
  particles <- cbind(c1, c2, kappa, sigma1, sigma2)
  
  print(apply(kappaStore[[iter]], 2, mean))
  print(apply(sigmaStore[[iter]], 2, mean))
  print(alpha)
}


