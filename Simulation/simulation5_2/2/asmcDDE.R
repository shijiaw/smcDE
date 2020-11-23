### adaptive annealing SMC, use B-spline to solve ODE###
ESSVec <- list()

parametersStore <- list()
sigmaStore <- list()
tempSeq <- list()
tempDifSeq <- list()
c1Store <- list()
c2Store <- list()
lambdaStore <- list()


logZ = 0
alpha <- 0
alphaDif <- 0

particles <- cbind(c1, c2, mu_m_particles, mu_p_particles, p_0_particles, tau_particles, sigma1, sigma2, lambda)
W <- rep(1/NP, NP)
logW <- log(W)
iter <- 0

while(alpha < 1){
  iter = iter + 1
  print(iter)
  c1 <- particles[,1:ncol(basismat)] 
  c2 <- particles[,(ncol(basismat)+1):(2*ncol(basismat))] 
  sigma1 <- particles[,2*ncol(basismat)+5] 
  sigma2 <- particles[,2*ncol(basismat)+6]
  sigma <- cbind(sigma1, sigma2)
  parameters <- particles[,(2*ncol(basismat)+1):(2*ncol(basismat)+4)]
  
  logRatio <- rep(NA, NP)
  for(i in 1:NP){
    logRatio[i] <- tempPosterior(as.vector(parameters[i,]), as.vector(sigma[i,]), y1, y2, as.vector(c1[i,]), as.vector(c2[i,]), lambda[i])
  }
  alphaDif <- bisection(func, 0, 1, W, logRatio, CESSthresholds)
  
  alpha <- alpha + alphaDif
  print(alpha)
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
  print(logZ)
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
  
  c1 <- particles[, 1:ncol(basismat)]
  c2 <- particles[, (ncol(basismat)+1):(2*ncol(basismat))]
  parameters <- particles[, (2*ncol(basismat)+1):(2*ncol(basismat)+4)]
  sigma1 <- particles[, (2*ncol(basismat)+5)]
  sigma2 <- particles[, (2*ncol(basismat)+6)]
  sigma <- cbind(sigma1, sigma2)
  for(i in 1:NP){
    ## Evaluate Basis functions ###
    Basis_eval <- newBasis(knots, parameters[i,4], bsbasis)
    quadpts_eval <- Basis_eval$quadpts
    quadwts_eval <- Basis_eval$quadwts
    D0quadbasismat_eval <- Basis_eval$D0quadbasismat
    D1quadbasismat_eval <- Basis_eval$D1quadbasismat
    D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
    
    ## Evaluate the curves ##
    xx1 <- basismat%*%c1[i,]
    xx2 <- basismat%*%c2[i,]
    dx1 <- D1quadbasismat_eval%*%c1[i,]
    dx2 <- D1quadbasismat_eval%*%c2[i,]
    x1 <- D0quadbasismat_eval%*%c1[i,]
    x2 <- D0quadbasismat_eval%*%c2[i,]
    tau_index = which(quadpts_eval == parameters[i,4])[1]
    x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2[i,])

    
    ##sample lambda using Gibbs
    lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, dx1, dx2, x1, x2,  x2_delay, quadwts_eval, alpha)
    
    ##sample sigma
    sigma[i,] <- GB_sigma(y1, y2, xx1, xx2, alpha)
    sigma1 <- sigma[,1]
    sigma2 <- sigma[,2]
    
    ##sample c1
    c1[i,] <- GB_c1(lambda[i], sigma1[i], parameters[i,], dx1, dx2, x1, x2,  x2_delay, quadwts_eval, D0quadbasismat_eval, D1quadbasismat_eval, y1, alpha)
    
    xx1 <- basismat%*%c1[i,]
    dx1 <- D1quadbasismat_eval%*%c1[i,]
    x1 <- D0quadbasismat_eval%*%c1[i,]
    
    ##sample c2
    c2_llOld <- c2_ll(xx2, y2, sigma2[i])
    c2_priorOld <- c2_logprior(lambda[i], parameters[i,], dx1, dx2, x1, x2, x2_delay, quadwts_eval, y2)
    MHstep_c2 <- MH_c2(c2_priorOld, c2_llOld, lambda[i], sigma2[i], parameters[i,], dx1, x1, c2[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y2, alpha)
    c2[i,] <- MHstep_c2$c2
    
    xx2 <- basismat%*%c2[i,]
    dx2 <- D1quadbasismat_eval%*%c2[i,]
    x2 <- D0quadbasismat_eval%*%c2[i,]
    x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2[i,])
    

    ##sample parameters
    parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i,], c1[i,], c2[i,], dx1, dx2, x1, x2, x2_delay, quadwts_eval, alpha)
    
  }
  proc.time() - ptm
  
  parametersStore[[iter]] <- parameters
  sigmaStore[[iter]] <- sigma
  c1Store[[iter]] <- c1
  c2Store[[iter]] <- c2
  lambdaStore[[iter]] <- lambda
  
  particles <- cbind(c1, c2, parameters, sigma1, sigma2, lambda)
  
  print(apply(parametersStore[[iter]], 2, mean))
  print(apply(sigmaStore[[iter]], 2, mean))
}





