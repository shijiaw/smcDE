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
    ##sample lambda using Gibbs
    lambda[i] <- GB_lambda(alambda, blambda, kappa[i,], M, dx1[,i], dx2[,i], x1[,i], x2[,i], alpha)
    
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
  
  #print(apply(kappaStore[[iter]], 2, mean))
  #print(apply(sigmaStore[[iter]], 2, mean))
  #print(alpha)
  #print(mean(lambda))
}


meankappa1 <- sum((kappaStore[[iter]][,1])*W)
kappa1Credible <- wtd.quantile(kappaStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
meankappa2 <- sum(kappaStore[[iter]][,2]*W)
kappa2Credible <- wtd.quantile(kappaStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
meansigma1 <- sum(sigmaStore[[iter]][,1]*W)
sigma1Credible <- wtd.quantile(sigmaStore[[iter]][,1], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))
meansigma2 <- sum(sigmaStore[[iter]][,2]*W)
sigma2Credible <- wtd.quantile(sigmaStore[[iter]][,2], na.rm = FALSE, weights = W, normwt=TRUE, prob = c(0.025, 0.975))



c1_estimate <- c1Store[[iter]]
c2_estimate <- c2Store[[iter]]

Spline_Matrix1 <- matrix(NA, nr = NP, nc = length(times))
Spline_Matrix2 <- matrix(NA, nr = NP, nc = length(times))
for(i in 1:NP){
  Spline_Matrix1[i,] <- basismat%*%c1_estimate[i,]
  Spline_Matrix2[i,] <- basismat%*%c2_estimate[i,]
}

MeanSpline1 <- Spline_Matrix1[1,]*W[1]
MeanSpline2 <- Spline_Matrix2[1,]*W[1] 
for(i in 2:NP){
  MeanSpline1 <- MeanSpline1 + Spline_Matrix1[i,]*W[i]
  MeanSpline2 <- MeanSpline2 + Spline_Matrix2[i,]*W[i]
}


c1Mean <- c1_estimate[1,]*W[1] 
c2Mean <- c2_estimate[1,]*W[1]

for(i in 2:nrow(c1_estimate)){ 
  c1Mean <- c1Mean + c1_estimate[i,]*W[i] 
  c2Mean <- c2Mean + c2_estimate[i,]*W[i] 
}


c1sd <- rep(0, ncol(c1_estimate)) #
c2sd <- rep(0, ncol(c2_estimate)) #
for(i in 1:nrow(c1_estimate)){ #   
  c1sd <- c1sd + (c1_estimate[i,] - c1Mean)^2*W[i] #   
  c2sd <- c2sd + (c2_estimate[i,] - c1Mean)^2*W[i] # 
} # #
c1sd <- sqrt(c1sd) # 
c2sd <- sqrt(c2sd) # # 
## plot ## # 
fitted1 <- basismat%*%c1Mean # 
upper1 <- fitted1 + 1.96*basismat%*%c1sd # 
lower1 <- fitted1 - 1.96*basismat%*%c1sd # # 
eval.basis(0, bsbasis)%*%c1Mean # 

fitted2 <- basismat%*%c2Mean # 
upper2 <- fitted2 + 1.96*basismat%*%c2sd # 
lower2 <- fitted2 - 1.96*basismat%*%c2sd # # 
eval.basis(0, bsbasis)%*%c2Mean #
plot(output[,2], type = 'l', ylim = c(-8, 11)) # 
lines(fitted1, col = 2) #
lines(upper1, lty = 3, col =2) # #
lines(lower1, lty = 3, col =2) # #
plot(output[,3], type = 'l', ylim = c(-80, 100)) # 
lines(fitted2, col = 2) #
lines(upper2, lty = 3, col =2) # #
lines(lower2, lty = 3, col =2) #

curve_matrix <- matrix(NA, nr = nrow(basismat), nc = NP)
for(i in 1:NP){
  curve_matrix[,i] <- basismat%*%c2_estimate[i,]
}
plot(curve_matrix[, 1], type = 'l', ylim = c(-30, 50))
for(i in 2:NP){
  lines(curve_matrix[, i], col = i)
}
upper_curve <- rep(NA, nrow(basismat))
mean_curve <- rep(NA, nrow(basismat))
lower_curve <- rep(NA, nrow(basismat))
for(j in 1:nrow(basismat)){
  upper_curve[j] <- wtd.quantile(curve_matrix[j,], probs = 1, type='quantile', weights=W, normwt=TRUE)
  lower_curve[j] <- wtd.quantile(curve_matrix[j,], probs = 0, type='quantile', weights=W, normwt=TRUE)
  mean_curve[j] <- wtd.quantile(curve_matrix[j,], probs = 0.5, type='quantile', weights=W, normwt=TRUE)
}
plot(mean_curve, type = 'l')
lines(output[,3], type = 'l', col = 3)
lines(upper_curve, type = 'l', col = 2)
lines(lower_curve, type = 'l', col = 2)
