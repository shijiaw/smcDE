###Metropolis Hastings algorithm ###
theta1_prop <- 4
theta2_prop <- 3

theta1_estimate <- rep(NA, nSamples)
theta2_estimate <- rep(NA, nSamples)

acc_rate <- 0

acc_event <- 0

MHsigma1 <- 0.1
MHsigma2 <- 0.1

theta1_estimate[1] <- 3
theta2_estimate[1] <- 2

l_old <- l(c(theta1_estimate[1], theta2_estimate[1]), y1, y2, sigma1, sigma2)

for(i in 2:nSamples){
  theta1_prop <- rnorm(1, theta1_estimate[i-1], MHsigma1)
  theta2_prop <- rnorm(1, theta2_estimate[i-1], MHsigma2)
  l_new <- l(c(theta1_prop, theta2_prop), y1, y2, sigma1, sigma2)
  log_prob <- l_new - l_old
  if(runif(1, 0, 1) <= exp(log_prob)){
    theta1_estimate[i] <- theta1_prop
    theta2_estimate[i] <- theta2_prop
    acc_event <- acc_event + 1
  }else{
    theta1_estimate[i] <- theta1_estimate[i-1]
    theta2_estimate[i] <- theta2_estimate[i-1]
  }
  print(i)
}




