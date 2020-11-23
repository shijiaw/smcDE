###Metropolis Hastings algorithm ###
thetaStore <- matrix(NA, nr = nSamples, nc = 2)
sigmaStore <- matrix(NA, nr = nSamples, nc = 2)
initialStore <- matrix(NA, nr = nSamples, nc = 2)

acc_rate <- 0

acc_event <- 0

theta <- c(4,3)
sigma <- c(2,5)
initial <- c(6, -11)

for(iter in 1:nSamples){
  print(iter)
  theta <- MHTheta(theta, initial, y1, y2, sigma)
  initial <- MHInitial(theta, initial, y1, y2, sigma)
  sigma <- GBKernel(theta, initial, y1, y2)
  thetaStore[iter,] <- theta
  sigmaStore[iter,] <- sigma
  initialStore[iter,] <- initial
}





