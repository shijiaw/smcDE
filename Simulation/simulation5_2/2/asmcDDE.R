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
  
  #if(iter > 1){
    #print(parametersStore[[iter]] - parametersStore[[iter-1]])
  #}
  
  particles <- cbind(c1, c2, parameters, sigma1, sigma2, lambda)
  
  print(apply(parametersStore[[iter]], 2, mean))
  print(apply(sigmaStore[[iter]], 2, mean))
  #print(mean(lambda[[iter]]))
}


meanmu_m <- sum(parametersStore[[iter]][,1]*W)
meanmu_p <- sum(parametersStore[[iter]][,2]*W)
meanp_0 <- sum(parametersStore[[iter]][,3]*W)
meantau <- sum(parametersStore[[iter]][,4]*W)
meansigma1 <- sum(sigmaStore[[iter]][,1]*W)
meansigma2 <- sum(sigmaStore[[iter]][,2]*W)
print(meanmu_m)
print(meanmu_p)
print(meanp_0)
print(meantau)
print(meansigma1)
print(meansigma2)


#sdmu_m <- sqrt(sum(W*(parametersStore[[iter]][,1] - meanmu_m)^2))
#sdmu_p <- sqrt(sum(W*(parametersStore[[iter]][,2] - meanmu_p)^2))
#print(c(meanmu_m-1.96*sdmu_m, meanmu_m+1.96*sdmu_m))
#print(c(meanmu_p-1.96*sdmu_p, meanmu_p+1.96*sdmu_p))

#sdp_0 <- sqrt(sum(W*(parametersStore[[iter]][,3] - meanp_0)^2))
#sdtau <- sqrt(sum(W*(parametersStore[[iter]][,4] - meantau)^2))
#print(c(meanp_0-1.96*sdp_0, meanp_0+1.96*sdp_0))
#print(c(meantau-1.96*sdtau, meantau+1.96*sdtau))


#sdsigma1 <- sqrt(sum(W*(sigmaStore[[iter]][,1] - meansigma1)^2))
#sdsigma2 <- sqrt(sum(W*(sigmaStore[[iter]][,2] - meansigma2)^2))
#print(c(meansigma1-1.96*sdsigma1, meansigma1+1.96*sdsigma1))
#print(c(meansigma2-1.96*sdsigma2, meansigma2+1.96*sdsigma2))

c1_estimate <- c1Store[[iter]]
c2_estimate <- c2Store[[iter]]
c1Mean <- c1_estimate[1,]*W[1]
c2Mean <- c2_estimate[1,]*W[1]

for(i in 2:nrow(c1_estimate)){
  c1Mean <- c1Mean + c1_estimate[i,]*W[i]
  c2Mean <- c2Mean + c2_estimate[i,]*W[i]
}

c1sd <- rep(0, ncol(c1_estimate))
c2sd <- rep(0, ncol(c2_estimate))
for(i in 1:nrow(c1_estimate)){
  c1sd <- c1sd + (c1_estimate[i,] - c1Mean)^2*W[i]
  c2sd <- c2sd + (c2_estimate[i,] - c2Mean)^2*W[i]
}

c1sd <- sqrt(c1sd)
c2sd <- sqrt(c2sd)

## plot ##
fitted1 <- basismat%*%c1Mean
upper1 <- fitted1 + 1.96*basismat%*%c1sd
lower1 <- fitted1 - 1.96*basismat%*%c1sd

eval.basis(0, bsbasis)%*%c1Mean

fitted2 <- basismat%*%c2Mean
upper2 <- fitted2 + 1.96*basismat%*%c2sd
lower2 <- fitted2 - 1.96*basismat%*%c2sd

eval.basis(0, bsbasis)%*%c2Mean
plot(output[,2], type = 'l', ylim = c(0, 30))
lines(fitted1, col = 2)
lines(upper1, lty = 3, col =2)
lines(lower1, lty = 3, col =2)

plot(output[,3], type = 'l', ylim = c(-80, 500))
lines(fitted2, col = 2)
lines(upper2, lty = 3, col =2)
lines(lower2, lty = 3, col =2)


#####code for processing results ###
mu_mCredible <- wtd.quantile(parametersStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(mu_mCredible)
mu_pCredible <- wtd.quantile(parametersStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(mu_pCredible)
p_0Credible <- wtd.quantile(parametersStore[[iter]][,3], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(p_0Credible)
tauCredible <- wtd.quantile(parametersStore[[iter]][,4], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(tauCredible)
sigma1Credible <- wtd.quantile(sigmaStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(sigma1Credible)
sigma2Credible <- wtd.quantile(sigmaStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(sigma2Credible)
C1Credible <- wtd.quantile(c1_estimate[,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(C1Credible)
C2Credible <- wtd.quantile(c2_estimate[,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
print(C2Credible)


### plot for DDE
c1_estimate <- c1Store[[iter]]
c2_estimate <- c2Store[[iter]]
c1Mean <- c1_estimate[1,]*W[1]
c2Mean <- c2_estimate[1,]*W[1]

for(i in 2:nrow(c1_estimate)){
  c1Mean <- c1Mean + c1_estimate[i,]*W[i]
  c2Mean <- c2Mean + c2_estimate[i,]*W[i]
}

curve_matrix1 <- matrix(NA, nr = nrow(basismat), nc = NP)
for(i in 1:NP){
  curve_matrix1[,i] <- basismat%*%c1_estimate[i,]
}

upper_curve1 <- rep(NA, nrow(basismat))
mean_curve1 <- rep(NA, nrow(basismat))
lower_curve1 <- rep(NA, nrow(basismat))
for(j in 1:nrow(basismat)){
  upper_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.975, type='quantile', weights=W, normwt=TRUE)
  lower_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.025, type='quantile', weights=W, normwt=TRUE)
  mean_curve1[j] <- wtd.quantile(curve_matrix1[j,], probs = 0.5, type='quantile', weights=W, normwt=TRUE)
}



curve_matrix2 <- matrix(NA, nr = nrow(basismat), nc = NP)
for(i in 1:NP){
  curve_matrix2[,i] <- basismat%*%c2_estimate[i,]
}

upper_curve2 <- rep(NA, nrow(basismat))
mean_curve2 <- rep(NA, nrow(basismat))
lower_curve2 <- rep(NA, nrow(basismat))
for(j in 1:nrow(basismat)){
  upper_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.975, type='quantile', weights=W, normwt=TRUE)
  lower_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.025, type='quantile', weights=W, normwt=TRUE)
  mean_curve2[j] <- wtd.quantile(curve_matrix2[j,], probs = 0.5, type='quantile', weights=W, normwt=TRUE)
}


plot(output[,3], type = 'l', ylim = c(-80, 500))
lines(mean_curve2, col = 2)
lines(upper_curve2, lty = 3, col =2)
lines(lower_curve2, lty = 3, col =2)

