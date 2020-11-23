##optimize ###
objectiveFunction4lambda <- function(par,c1,c2){
  kappa1 <- par[1]
  kappa2 <- par[2]
  rt1 <- sum((y1 - basismat%*%c1)^2/1^2+(y2 - basismat%*%c2)^2/2^2)
  dx1 <- D1quadbasismat%*%c1
  dx2 <- D1quadbasismat%*%c2
  x1 <- D0quadbasismat%*%c1
  x2 <- D0quadbasismat%*%c2
  g1 <- 72/(36+x2) - kappa1 
  g2 <- kappa2*x1 - 1
  rt2 <- sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
  return(rt2)
}

 
par_parameter <- c(5, 7)
optimizer4lambda <- nlminb(start = par_parameter, objectiveFunction4lambda, control=list(trace=FALSE, abs.tol = 0.1, iter.max = 100), c1 = c1_LSE, c2 = c2_LSE)
initialized <- optimizer4lambda$par


dx1 <- D1quadbasismat%*%c1_LSE
dx2 <- D1quadbasismat%*%c2_LSE
x1 <- D0quadbasismat%*%c1_LSE
x2 <- D0quadbasismat%*%c2_LSE
g1 <- 72/(36+x2) - initialized[1] 
g2 <- initialized[2]*x1 - 1
lambda1 <- sum(quadwts*((dx1-g1)^2))
lambda2 <- sum(quadwts*((dx2-g2)^2))
#rt2 <- sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
#lambda <- sqrt(min(lambda1, lambda2))

lambda <- 1
