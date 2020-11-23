### adaptive annealing SMC, use B-spline to solve ODE###
ESSVec <- list()

parametersStore <- list()
sigmaStore <- list()
#initialStore <- list()
tempSeq <- list()
tempDifSeq <- list()
cStore <- list()
lambdaStore <- list()


logZ = 0
alpha <- 0
alphaDif <- 0
#tempDif <- tempScheme[-1] - tempScheme[-T]
#tempDif <- c(tempScheme[1], tempDif)

#particles <- cbind(theta1, theta2, sigma1, sigma2, c1, c2, lambda)
#particles <- list(theta1 = theta1, theta2 = theta2, sigma1 = sigma1, sigma2 = sigma2, c1 = c1, c2 = c2, lambda = lambda)
particles <- cbind(c, r_particles, P_particles, tau_particles, sigma, lambda)
### compute the first temperature ###
#logw <- apply(as.matrix(particles), 1, function(x)(Loglikelihood(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))
#logRatio <- apply(as.matrix(particles), 1, function(x)(Loglikelihood(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))


### double check weights for normalized version ###
W <- rep(1/NP, NP)
logW <- log(W)
iter <- 0
#logW <- log(W)

#logmax = max(logw)
#W = exp(logw-logmax)/sum(exp(logw-logmax)) 
#RESS <- ESS(W)

while(alpha < 1){
  iter = iter + 1
  print(iter)
  ### Compute weights and conduct resampling ###
  #logRatio <- apply(as.matrix(particles), 1, function(x)(Loglikelihood(as.vector(x)[1:2], as.vector(x)[3:4], as.vector(x)[5:6], y1, y2)))
  c <- particles[,1:ncol(basismat)] 
  sigma <- particles[,ncol(basismat)+4] 
  parameters <- particles[,(ncol(basismat)+1):(ncol(basismat)+3)]
  #lambda <- rep(1, NP)
  logRatio <- rep(NA, NP)
  for(i in 1:NP){
    logRatio[i] <- tempPosterior(as.vector(parameters[i,]), sigma[i], y, as.vector(c[i,]), lambda[i])
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
    #ancestor = sample.int(NP, prob = W, replace = TRUE)
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
    ## Evaluate Basis functions ###
    Basis_eval <- newBasis(knots, parameters[i,3], bsbasis)
    quadpts_eval <- Basis_eval$quadpts
    quadwts_eval <- Basis_eval$quadwts
    D0quadbasismat_eval <- Basis_eval$D0quadbasismat
    D1quadbasismat_eval <- Basis_eval$D1quadbasismat
    D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
    
    ## Evaluate the curves ##
    ww <- basismat%*%c[i,]
    dw <- D1quadbasismat_eval%*%c[i,]
    w <- D0quadbasismat_eval%*%c[i,]
    tau_index = which(quadpts_eval == parameters[i,3])[1]
    w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
    
    
    ##sample lambda using Gibbs
    lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, w, dw, w_delay, quadwts_eval, alpha)
    
    ##sample sigma
    sigma[i] <- GB_sigma(y, ww, alpha)
    
    ##sample c
    c_llOld <- c_ll(ww, y, sigma[i])
    c_priorOld <- c_logprior(lambda[i], parameters[i,], dw, w_delay, quadwts_eval)
    MHstep_c <- MH_c(c_priorOld, c_llOld, lambda[i], sigma[i], parameters[i,], c[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y, alpha)
    c[i,] <- MHstep_c$c
    
    ww <- basismat%*%c[i,]
    dw <- D1quadbasismat_eval%*%c[i,]
    w <- D0quadbasismat_eval%*%c[i,]
    tau_index = which(quadpts_eval == parameters[i,3])[1]
    w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
    
    
    ##sample parameters
    parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i], c[i,], dw, w, w_delay, quadwts_eval, alpha)

    
  }
  proc.time() - ptm
  
  
  parametersStore[[iter]] <- parameters
  sigmaStore[[iter]] <- sigma
  cStore[[iter]] <- c
  lambdaStore[[iter]] <- lambda
  
  if(iter > 1){
    #print(parametersStore[[iter]] - parametersStore[[iter-1]])
  }
  
  particles <- cbind(c, parameters, sigma, lambda)
  
  print(apply(parametersStore[[iter]], 2, mean))
  print(mean(sigmaStore[[iter]]))
  #print(mean(lambda[[iter]]))
}

mean_r <- sum(parametersStore[[iter]][,1]*W)
mean_P <- sum(parametersStore[[iter]][,2]*W)
mean_tau <- sum(parametersStore[[iter]][,3]*W)
meansigma <- sum(sigmaStore[[iter]]^2*W)

mean_Q <- sum(parametersStore[[iter]][,1]/parametersStore[[iter]][,2]*W)

### plot for DDE
c_estimate <- cStore[[iter]]

cMean <- c_estimate[1,]*W[1]


for(i in 2:nrow(c_estimate)){
  cMean <- cMean + c_estimate[i,]*W[i]
}

csd <- rep(0, ncol(c_estimate))

for(i in 1:nrow(c_estimate)){
  csd <- csd + (c_estimate[i,] - cMean)^2*W[i]
}

csd <- sqrt(csd)


## plot ##
fitted1 <- basismat%*%cMean
upper1 <- fitted1 + 1.96*basismat%*%csd
lower1 <- fitted1 - 1.96*basismat%*%csd

eval.basis(0, bsbasis)%*%cMean



plot(log(output[,2]), type = 'l', ylim = c(0, 20))
lines((fitted1), col = 2)
lines((upper1), lty = 3, col =2)
lines((lower1), lty = 3, col =2)


#####code for processing results ###
r_Credible <- wtd.quantile(parametersStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(r_Credible)
P_Credible <- wtd.quantile(parametersStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(P_Credible)
tau_Credible <- wtd.quantile(parametersStore[[iter]][,3], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(tau_Credible)
sigma_Credible <- wtd.quantile(sigmaStore[[iter]]^2, na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(sigma_Credible)

cMean[1]
W_Credible <- wtd.quantile(c_estimate[,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
sqrt(sum((log(output[,2]) - fitted1)^2)/length(fitted1))
wtd.quantile(sigmaStore[[iter]]^2, na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))

