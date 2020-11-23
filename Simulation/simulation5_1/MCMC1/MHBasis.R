### MCMC, use B-spline to solve ODE###
ptm <- proc.time()
for(iter in 2:Niter){

  lambda <- GB_lambda(alambda, blambda, kappa, M, c1, c2)
  
  ##sample sigma
  sigma <- GB_sigma(y1, y2, c1, c2)
  
  ##sample kappa
  kappa <- MH_kappa(lambda, kappa, sigma, c1, c2)
  
  ##sample c1
  c1 <- GB_c1(lambda, sigma[1], kappa, c1, c2, y1)
  
  ##sample c2
  c2_llOld <- c2_ll(c2, y2, sigma[2])
  c2_priorOld <- c2_logprior(lambda, kappa, c1, c2, y2)
  MHstep_c2 <- MH_c2(c2_priorOld, c2_llOld, lambda, sigma[2], kappa, c1, c2, y2)
  c2 <- MHstep_c2$c2
  
  kappaStore[iter,] <- kappa
  sigmaStore[iter,] <- sigma
  c1Store[iter,] <- c1
  c2Store[iter,] <- c2
  lambdaStore[iter] <- lambda
  print(iter)
}
proc.time() - ptm
RMSE1 <- sqrt(sum((basismat%*%apply(c1Store[-(1:200000),], 2, mean) - true1)^2)/length(true1))
RMSE2 <- sqrt(sum((basismat%*%apply(c2Store[-(1:200000),], 2, mean) - true2)^2)/length(true2))


filename_theta1 = './theta1'
filename_theta2 = './theta2'
filename_sigma1 = './sigma1'
filename_sigma2 = './sigma2'
filename_lambda = './lambda'


write.table(matrix(round(kappaStore[,1], digits = 5),nr=1), file = paste(filename_theta1,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(kappaStore[,2], digits = 5),nr=1), file = paste(filename_theta2,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(sigmaStore[,1], digits = 5),nr=1), file = paste(filename_sigma1,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(sigmaStore[,2], digits = 5),nr=1), file = paste(filename_sigma2,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(lambdaStore, digits = 5),nr=1), file = paste(filename_lambda,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)


