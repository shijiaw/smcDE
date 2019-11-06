#' This function estimate parameters in DEs via SMC based on collocation 
#' @param data a list with each element being a trajectory with measurement error
#' @param times a list with each element being the time we record data
#' @param seed the random seeds for propagation and resampling in SMC
#' @param knotsPosition the location we put knots
#' @param NP the number of particles in SMC
#' @param CESSthresholds relative conditional effective sample size (0, 1) in SMC
#' @param resampleThreshold the threshold triggering resampling in SMC
#' @param alambda hyper-parameter for tuning parameter lambda
#' @param blambda hyper-parameter for tuning parameter lambda
#' @param sigmac the standard deviation of the basis coefficients c
#' @param DEmodel options for DE models: 1 (ODE), 2 (DDE1), 3 (DDE2) in manuscript
#' @export 

smcDE <- function(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac, DEmodel){
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  # check if packages are installed
  packages<-c("fda", "truncnorm", "numDeriv", "deSolve", "NLRoot", "smcUtils","Hmisc","MASS")
  check.packages(packages)
  
  # load required packages
  library(truncnorm)
  library(numDeriv)
  library(fda)
  library(deSolve)
  library(NLRoot)
  library(smcUtils)
  library(Hmisc)
  library(MASS)
  
  # Three DE models: ODE, DDE1, DDE2
  if((DEmodel != 1) && (DEmodel != 2) && (DEmodel != 3)){
    cat("The input (DEmodel) is invalid! Three options: 1 (ODE), 2 (DDE1), 3 (DDE2)", "\n")
  }
  
  if(length(data) != length(times)){
    cat("Error in the input data and times", "\n")
  }
  
  dir <- getwd()
  
  if(DEmodel == 1){
    #source("DE1.R")
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
    
    output <- DE1(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  if(DEmodel == 2){
    #source("DE2.R")
    DE2 <- function(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac){
      y = data[[1]]
      times = times[[1]]
      
      # set random seed for smc 
      set.seed(seed)
      
      tobs = times
      knots = knots
      nknots   = length(knots)
      norder   = 4
      nbasis   = length(knots) + norder - 2
      bsbasis = create.bspline.basis(range(knots), nbasis, norder, knots)
      
      # basis values at sampling points
      basismat = eval.basis(tobs, bsbasis)
      
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
      
      M <- nknots
      
      # initialize particles 
      # Note that the particles can also be generated from other distributions (e.g. priors)
      r_particles = runif(NP, 0, 2)
      P_particles = runif(NP, 0, 5)
      tau_particles = runif(NP, 0, 10)
      
      sigma <- runif(NP, 0.3, 5)
      
      lambda <- rgamma(NP, 1, 10)
      
      c_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%log(y)
      
      c <- matrix(NA, nr = NP, nc = ncol(basismat))
      
      for(i in 1:NP){
        c[i,] <- rnorm(length(c_LSE), c_LSE, 0.2)
      }
      
      # adaptive annealing SMC, use B-spline to solve ODE
      ESSVec <- list()
      parametersStore <- list()
      sigmaStore <- list()
      tempSeq <- list()
      tempDifSeq <- list()
      cStore <- list()
      lambdaStore <- list()
      
      
      logZ = 0
      alpha <- 0
      alphaDif <- 0
      
      particles <- cbind(c, r_particles, P_particles, tau_particles, sigma, lambda)
      
      # weights function 
      W <- rep(1/NP, NP)
      logW <- log(W)
      iter <- 0
      
      while(alpha < 1){
        iter = iter + 1
        cat("iteration:", iter, "\n")
        
        # Compute weights and conduct resampling 
        c <- particles[,1:ncol(basismat)] 
        sigma <- particles[,ncol(basismat)+4] 
        parameters <- particles[,(ncol(basismat)+1):(ncol(basismat)+3)]
        
        logRatio <- rep(NA, NP)
        
        for(i in 1:NP){
          logRatio[i] <- tempPosterior(as.vector(parameters[i,]), sigma[i], y, as.vector(c[i,]), lambda[i], basismat, knots, bsbasis, nquad)
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
        
        c <- particles[, 1:ncol(basismat)]
        parameters <- particles[, (ncol(basismat)+1):(ncol(basismat)+3)]
        sigma <- particles[, (ncol(basismat)+4)]
        
        for(i in 1:NP){
          # Evaluate Basis functions 
          Basis_eval <- newBasis(knots, parameters[i,3], bsbasis, nquad)
          quadpts_eval <- Basis_eval$quadpts
          quadwts_eval <- Basis_eval$quadwts
          D0quadbasismat_eval <- Basis_eval$D0quadbasismat
          D1quadbasismat_eval <- Basis_eval$D1quadbasismat
          D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
          
          # Evaluate the curves 
          ww <- basismat%*%c[i,]
          dw <- D1quadbasismat_eval%*%c[i,]
          w <- D0quadbasismat_eval%*%c[i,]
          tau_index = which(quadpts_eval == parameters[i,3])[1]
          w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
          
          # sample lambda using Gibbs
          lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, w, dw, w_delay, quadwts_eval, alpha)
          
          # sample sigma
          sigma[i] <- GB_sigma(y, ww, alpha)
          
          # sample c
          c_llOld <- c_ll(ww, y, sigma[i])
          c_priorOld <- c_logprior(lambda[i], parameters[i,], dw, w_delay, quadwts_eval)
          MHstep_c <- MH_c(c_priorOld, c_llOld, lambda[i], sigma[i], parameters[i,], c[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y, alpha)
          c[i,] <- MHstep_c$c
          
          ww <- basismat%*%c[i,]
          dw <- D1quadbasismat_eval%*%c[i,]
          w <- D0quadbasismat_eval%*%c[i,]
          tau_index = which(quadpts_eval == parameters[i,3])[1]
          w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c[i,])
          
          # sample parameters
          parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i], c[i,], dw, w, w_delay, quadwts_eval, alpha, knots, bsbasis, nquad)
        }
        
        proc.time() - ptm
        
        parametersStore[[iter]] <- parameters
        sigmaStore[[iter]] <- sigma
        cStore[[iter]] <- c
        lambdaStore[[iter]] <- lambda
        
        particles <- cbind(c, parameters, sigma, lambda)
      }
      
      # particles at final iteration
      c <- particles[, 1:ncol(basismat)]
      parameters <- particles[, (ncol(basismat)+1):(ncol(basismat)+3)]
      sigma <- particles[, (ncol(basismat)+4)]
      lambda <- particles[, (ncol(basismat)+5)]
      
      parameterList <- list(r = parameters[,1], P = parameters[,2], tau = parameters[,3])
      
      
      return(list(c = c, parameters = parameterList, lambda = lambda, sigma = sigma, W = W))
    }
    
    # function for evaluating basis functions with new tau sampled
    newBasis <- function(knots, tau, bsbasis, nquad){
      new_eval = sort(c(knots, tau))
      quadre = quadset(nquad, bsbasis, new_eval)
      quadpts = quadre$quadvals[,1]
      
      tau_index = which(quadpts == tau)[1]
      
      D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
      D0quadbasismat_delay = eval.basis((quadpts-tau)[-(1:tau_index)], bsbasis, 0)
      D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
      return(list(quadpts=quadpts, quadwts=quadre$quadvals[,2], D0quadbasismat = D0quadbasismat, D1quadbasismat = D1quadbasismat, D0quadbasismat_delay = D0quadbasismat_delay))
    }
    
    # Compute tempering posterior
    tempPosterior <- function(parameters, sigma, y, c, lambda, basismat, knots, bsbasis, nquad){
      y <- log(y)
      
      r <- parameters[1]
      P <- parameters[2]
      tau <- parameters[3]
      
      sigma <- as.numeric(sigma)
      
      Basis_eval <- newBasis(knots, tau, bsbasis, nquad)
      quadpts_eval <- Basis_eval$quadpts
      quadwts_eval <- Basis_eval$quadwts
      D0quadbasismat_eval <- Basis_eval$D0quadbasismat
      D1quadbasismat_eval <- Basis_eval$D1quadbasismat
      D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
      
      ww <- basismat%*%c
      
      dw <- D1quadbasismat_eval%*%c
      w <- D0quadbasismat_eval%*%c
      
      tau_index = which(quadpts_eval == tau)[1]
      w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%c)
      
      g <- r*(1-exp(w_delay)/(1000*P))
      
      rt <- -sum((y - ww)^2/sigma^2) -length(y)*log(sigma^2) - lambda*sum(quadwts_eval*((dw-g)^2))
      return(rt)
    }
    
    
    # Metropolis Hastings sampling for parameters (r, p, tau)
    MH_parameter <- function(lambda, parametersOld, sigma, c, dw, w, w_delay, quadwts, temperature, knots, bsbasis, nquad){
      r_Old <- parametersOld[1]
      P_Old <- parametersOld[2]
      tauOld <- parametersOld[3]
      Q_Old <- r_Old/P_Old
      
      sigmar_Inv = temperature*lambda*sum(quadwts)
      sigmar = 1/sqrt(sigmar_Inv)
      mu_r = temperature*lambda*(sum(quadwts*(dw+Q_Old*exp(w_delay)/1000)))/sigmar_Inv
      rnew <- rtruncnorm(1, a=0, b=Inf, mean = mu_r, sd = sigmar)
      
      sigmaQ_Inv <- temperature*lambda*sum(quadwts*(exp(w_delay)/(1000))^2)
      sigmaQ <- 1/sqrt(sigmaQ_Inv)
      mu_Q <- -temperature*lambda*(sum(quadwts*(dw-rnew)*(exp(w_delay)/(1000))))/sigmaQ_Inv
      Qnew <- rtruncnorm(1, a=0, b=Inf, mean = mu_Q, sd = sigmaQ)
      P_new <- rnew/Qnew
      
      g1Old <- rnew*(1-exp(w_delay)/(1000*P_new))
      
      u2 <- runif(1)
      if(u2 < 0.5){
        tausigma <- 0.2
      }
      if(u2 >= 0.5){
        tausigma <- 2
      }
      # taunew <- rnorm(1, tauOld, tausigma)
      taunew <- rtruncnorm(1, a=0, b=10, mean = tauOld, sd = tausigma)
      
      ## re-define basis eval ##
      Basis_eval <- newBasis(knots, taunew, bsbasis, nquad)
      quadpts_eval <- Basis_eval$quadpts
      quadwts_eval <- Basis_eval$quadwts
      D0quadbasismat_eval <- Basis_eval$D0quadbasismat
      D1quadbasismat_eval <- Basis_eval$D1quadbasismat
      D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
      
      dw_new <- D1quadbasismat_eval%*%c
      w_new <- D0quadbasismat_eval%*%c
      tau_index = which(quadpts_eval == taunew)[1]
      w_delay <- c(rep(w_new[1], tau_index), D0quadbasismat_delay%*%c)
      
      g1new <- rnew*(1-exp(w_delay)/(1000*P_new))
      
      logacc_prob <- temperature*(-lambda*sum(quadwts_eval*((dw_new-g1new)^2)) + lambda*sum(quadwts*((dw-g1Old)^2)))
      if(runif(1, 0, 1) < exp(logacc_prob)){
        rt <- c(rnew, P_new, taunew)
      }else{
        rt <- c(rnew, P_new, tauOld)
      }  
      
      return(rt)
    }
    
    # log likelihood for c
    c_ll <- function(ww, y, sigma){
      y <- log(y)
      #w <- basismat%*%c
      rt <- -sum((y - ww)^2/sigma^2)
      return(rt)
    }
    
    # prior distribution for c
    c_logprior <- function(lambda, parameters, dw, w_delay, quadwts){
      r <- parameters[1]
      P <- parameters[2]
      tau <- parameters[3]
      
      g <- r*(1-exp(w_delay)/(1000*P))
      
      rt <- -lambda*sum(quadwts*((dw-g)^2))
      return(rt)
    }
    
    # Metropolis Hasting to sample basis coefficients c
    MH_c <- function(c_priorOld, c_llOld, lambda, sigma, parameters, c, quadwts, tau_index, basismat, D0quadbasismat, D1quadbasismat, D0quadbasismat_delay, y, temperature){
      #MHsigma <- 1
      newc <- c(rnorm(length(c), c, sigmac))
      
      ww <- basismat%*%newc
      dw <- D1quadbasismat%*%newc
      w <- D0quadbasismat%*%newc
      w_delay <- c(rep(w[1], tau_index), D0quadbasismat_delay%*%newc)
      
      c_priornew <- c_logprior(lambda, parameters, dw, w_delay, quadwts)
      c_llnew <- c_ll(ww, y, sigma)
      logacc_prob <- temperature*(c_priornew - c_priorOld) + temperature*(c_llnew - c_llOld)
      if(runif(1, 0, 1) < exp(logacc_prob)){
        #print("accept")
        rt <- c(newc)
      }else{
        #print("reject")
        rt <- c(c)
      }
      return(list(c_prior = c_priornew, c_ll = c_llnew, c =rt)) 
    }
    
    # Gibbs step for measurement error sigma 
    GB_sigma <- function(y, ww, temperature){
      y <- log(y)
      temp <- rgamma(1, shape = 1+ length(y)/2, rate = 1 + temperature*sum((y - ww)^2)/2)
      return(sqrt(1/temp))
    }
    
    # Gibbs move for lambda 
    GB_lambda <- function(alambda, blambda, parameters, M, w, dw,  w_delay, quadwts, temperature){
      r <- parameters[1]
      P <- parameters[2]
      tau <- parameters[3]
      g <- r*(1-exp(w_delay)/(1000*P))
      b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dw-g)^2))/2)
      return(rgamma(1, shape = alambda+temperature*M/2, scale = b_lambda))
    }
    
    
    # effective sample size
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
    
    
    output <- DE2(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  if(DEmodel == 3){
    #source(("DE3.R"))
    DE3 <- function(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac){
      y1 = data[[1]]
      y2 = data[[2]]
      times = times[[1]]
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
      
      M <- nknots
      
      mu_m_particles = abs(rnorm(NP, 0.04, 0.02))
      mu_p_particles = abs(rnorm(NP, 0.04, 0.02))
      p_0_particles = rnorm(NP, 90, 20)
      tau_particles = rnorm(NP, 30, 10)
      
      sigma1 <- runif(NP, 0.3, 3)
      sigma2 <- runif(NP, 0.2, 8)
      
      lambda <- (rgamma(NP, 1, 1))
      
      c1_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y1
      c2_LSE <- solve(as.matrix(basismat2))%*%t(basismat)%*%y2
      
      double_quadpts <- table(quadpts)[which(table(quadpts) == 2)]
      
      c1 <- matrix(NA, nr = NP, nc = ncol(basismat))
      c2 <- matrix(NA, nr = NP, nc = ncol(basismat))
      
      for(i in 1:NP){
        c1[i,] <- rnorm(length(c1_LSE), c1_LSE, 0.1)
        c2[i,] <- rnorm(length(c2_LSE), c2_LSE, 0.1)
      }
      
      # adaptive annealing SMC, use B-spline to solve ODE 
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
      
      # weights function 
      W <- rep(1/NP, NP)
      logW <- log(W)
      iter <- 0
      
      while(alpha < 1){
        iter = iter + 1
        cat("iteration:", iter, "\n")
        # Compute weights and conduct resampling 
        c1 <- particles[,1:ncol(basismat)] 
        c2 <- particles[,(ncol(basismat)+1):(2*ncol(basismat))] 
        sigma1 <- particles[,2*ncol(basismat)+5] 
        sigma2 <- particles[,2*ncol(basismat)+6]
        sigma <- cbind(sigma1, sigma2)
        parameters <- particles[,(2*ncol(basismat)+1):(2*ncol(basismat)+4)]
        
        logRatio <- rep(NA, NP)
        for(i in 1:NP){
          logRatio[i] <- tempPosterior(as.vector(parameters[i,]), as.vector(sigma[i,]), y1, y2, as.vector(c1[i,]), as.vector(c2[i,]), lambda[i], basismat, knots, bsbasis, nquad)
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
        
        c1 <- particles[, 1:ncol(basismat)]
        c2 <- particles[, (ncol(basismat)+1):(2*ncol(basismat))]
        parameters <- particles[, (2*ncol(basismat)+1):(2*ncol(basismat)+4)]
        sigma1 <- particles[, (2*ncol(basismat)+5)]
        sigma2 <- particles[, (2*ncol(basismat)+6)]
        sigma <- cbind(sigma1, sigma2)
        for(i in 1:NP){
          # Evaluate Basis functions 
          Basis_eval <- newBasis(knots, parameters[i,4], bsbasis, nquad)
          quadpts_eval <- Basis_eval$quadpts
          quadwts_eval <- Basis_eval$quadwts
          D0quadbasismat_eval <- Basis_eval$D0quadbasismat
          D1quadbasismat_eval <- Basis_eval$D1quadbasismat
          D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
          
          # Evaluate the curves 
          xx1 <- basismat%*%c1[i,]
          xx2 <- basismat%*%c2[i,]
          dx1 <- D1quadbasismat_eval%*%c1[i,]
          dx2 <- D1quadbasismat_eval%*%c2[i,]
          x1 <- D0quadbasismat_eval%*%c1[i,]
          x2 <- D0quadbasismat_eval%*%c2[i,]
          tau_index = which(quadpts_eval == parameters[i,4])[1]
          x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2[i,])
          
          
          # sample lambda using Gibbs
          lambda[i] <- GB_lambda(alambda, blambda, parameters[i,], M, dx1, dx2, x1, x2,  x2_delay, quadwts_eval, alpha)
          
          # sample sigma
          sigma[i,] <- GB_sigma(y1, y2, xx1, xx2, alpha)
          sigma1 <- sigma[,1]
          sigma2 <- sigma[,2]
          
          # sample c1
          c1[i,] <- GB_c1(lambda[i], sigma1[i], parameters[i,], dx1, dx2, x1, x2,  x2_delay, quadwts_eval, D0quadbasismat_eval, D1quadbasismat_eval, y1, alpha, basismat)
          
          xx1 <- basismat%*%c1[i,]
          dx1 <- D1quadbasismat_eval%*%c1[i,]
          x1 <- D0quadbasismat_eval%*%c1[i,]
          
          # sample c2
          c2_llOld <- c2_ll(xx2, y2, sigma2[i])
          c2_priorOld <- c2_logprior(lambda[i], parameters[i,], dx1, dx2, x1, x2, x2_delay, quadwts_eval, y2)
          MHstep_c2 <- MH_c2(c2_priorOld, c2_llOld, lambda[i], sigma2[i], parameters[i,], dx1, x1, c2[i,], quadwts_eval, tau_index, basismat, D0quadbasismat_eval, D1quadbasismat_eval, D0quadbasismat_delay, y2, alpha, sigmac)
          c2[i,] <- MHstep_c2$c2
          
          xx2 <- basismat%*%c2[i,]
          dx2 <- D1quadbasismat_eval%*%c2[i,]
          x2 <- D0quadbasismat_eval%*%c2[i,]
          x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2[i,])
          
          # sample parameters
          parameters[i,] <- MH_parameter(lambda[i], parameters[i,], sigma[i,], c1[i,], c2[i,], dx1, dx2, x1, x2, x2_delay, quadwts_eval, alpha, knots, bsbasis, nquad)
          
        }
        proc.time() - ptm
        
        parametersStore[[iter]] <- parameters
        sigmaStore[[iter]] <- sigma
        c1Store[[iter]] <- c1
        c2Store[[iter]] <- c2
        lambdaStore[[iter]] <- lambda
        
        particles <- cbind(c1, c2, parameters, sigma1, sigma2, lambda)
      }
      
      parameterList <- list(mu_m = parameters[,1], mu_p = parameters[,2], p_0 = parameters[,3], tau = parameters[,4])
      
      return(list(c1 = c1, c2 = c2, parameters = parameterList, lambda = lambda, sigma = sigma, W = W))
    }
    
    # function for evaluating basis functions with new tau sampled
    newBasis <- function(knots, tau, bsbasis, nquad){
      new_eval = sort(c(knots, tau))
      quadre = quadset(nquad, bsbasis, new_eval)
      quadpts=quadre$quadvals[,1]
      
      tau_index = which(quadpts == tau)[1]
      
      D0quadbasismat = eval.basis(quadpts, bsbasis, 0)
      D0quadbasismat_delay = eval.basis((quadpts-tau)[-(1:tau_index)], bsbasis, 0)
      D1quadbasismat = eval.basis(quadpts, bsbasis, 1)
      return(list(quadpts=quadpts, quadwts=quadre$quadvals[,2], D0quadbasismat = D0quadbasismat, D1quadbasismat = D1quadbasismat, D0quadbasismat_delay = D0quadbasismat_delay))
    }
    
    
    tempPosterior <- function(parameters, sigma, y1, y2, c1, c2, lambda, basismat, knots, bsbasis, nquad){
      mu_m <- parameters[1]
      mu_p <- parameters[2]
      p_0 <- parameters[3]
      tau <- parameters[4]
      
      sigma1 <- as.numeric(sigma[1])
      sigma2 <- as.numeric(sigma[2])
      
      Basis_eval <- newBasis(knots, tau, bsbasis, nquad)
      quadpts_eval <- Basis_eval$quadpts
      quadwts_eval <- Basis_eval$quadwts
      D0quadbasismat_eval <- Basis_eval$D0quadbasismat
      D1quadbasismat_eval <- Basis_eval$D1quadbasismat
      D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
      
      xx1 <- basismat%*%c1
      xx2 <- basismat%*%c2
      dx1 <- D1quadbasismat_eval%*%c1
      dx2 <- D1quadbasismat_eval%*%c2
      x1 <- D0quadbasismat_eval%*%c1
      x2 <- D0quadbasismat_eval%*%c2
      tau_index = which(quadpts_eval == tau)[1]
      x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%c2)
      
      g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
      g2 <- x1 - mu_p*x2
      
      rt <- -sum((y1 - xx1)^2/sigma1^2 + (y2 - xx2)^2/sigma2^2) -length(y1)*log(sigma1^2) -length(y2)*log(sigma2^2) - lambda*sum(quadwts_eval*((dx1-g1)^2+(dx2-g2)^2))
      return(rt)
    }
    
    
    # Metropolis Hastings sampling for parameters mu_p, mu_p, p_0, tau
    MH_parameter <- function(lambda, parametersOld, sigma, c1, c2, dx1, dx2, x1, x2, x2_delay, quadwts, temperature, knots, bsbasis, nquad){
      mu_mOld <- parametersOld[1]
      mu_pOld <- parametersOld[2]
      p_0Old <- parametersOld[3]
      tauOld <- parametersOld[4]
      
      temp_delay <- (x2_delay/p_0Old)^DDE_power
      
      sigmamu_mInv = temperature*lambda*sum(quadwts*x1^2)
      sigmamu_m = 1/sqrt(sigmamu_mInv)
      mu_mu_m = -temperature*lambda*(sum(quadwts*x1*(dx1-1/(1+temp_delay))))/sigmamu_mInv
      mu_mnew <- rtruncnorm(1, a=0, b=100, mean = mu_mu_m, sd = sigmamu_m)
      
      sigmamu_pInv = temperature*lambda*sum(quadwts*x2^2)
      sigmamu_p = 1/sqrt(sigmamu_pInv)
      mu_mu_p = temperature*lambda*(sum(quadwts*x2*(x1-dx2)))/sigmamu_pInv
      mu_pnew <- rtruncnorm(1, a=0, b=100, mean = mu_mu_p, sd = sigmamu_p)
      
      u1 <- runif(1)
      if(u1 < 0.5){
        psigma <- 2
      }
      if(u1 >= 0.5){
        psigma <- 10
      }
      
      p_0new <- rtruncnorm(1, a=0, b=Inf, mean = p_0Old, sd = psigma)
      
      g1Old <- 1/(1+ (x2_delay/p_0Old)^DDE_power) - mu_mnew*x1 
      g1new <- 1/(1+ (x2_delay/p_0new)^DDE_power) - mu_mnew*x1 
      #g2new <- x1 - mu_pnew*x2
      
      logacc_prob <- dnorm(p_0new, 90, 20, log = TRUE) - dnorm(p_0Old, 90, 20, log = TRUE) +temperature*(-lambda*sum(quadwts*((dx1-g1new)^2)) + lambda*sum(quadwts*((dx1-g1Old)^2)))
      #logacc_prob <- tempLfun(newTheta, newInitial, sigma, y1, y2, temperature) - tempLfun(oldTheta, oldInitial, sigma, y1, y2, temperature)
      if(runif(1, 0, 1) > exp(logacc_prob)){
        p_0new <- p_0Old
      }
      
      g1Old <- 1/(1+(x2_delay/p_0new)^DDE_power) - mu_mnew*x1 
      
      u2 <- runif(1)
      if(u2 < 0.5){
        tausigma <- 1
      }
      if(u2 >= 0.5){
        tausigma <- 5
      }
      taunew <- rtruncnorm(1, a=0, b=50, mean = tauOld, sd = tausigma)
      
      ## re-define basis eval ##
      Basis_eval <- newBasis(knots, taunew, bsbasis, nquad)
      quadpts_eval <- Basis_eval$quadpts
      quadwts_eval <- Basis_eval$quadwts
      D0quadbasismat_eval <- Basis_eval$D0quadbasismat
      D1quadbasismat_eval <- Basis_eval$D1quadbasismat
      D0quadbasismat_delay <- Basis_eval$D0quadbasismat_delay
      
      dx1_new <- D1quadbasismat_eval%*%c1
      dx2_new <- D1quadbasismat_eval%*%c2
      x1_new <- D0quadbasismat_eval%*%c1
      x2_new <- D0quadbasismat_eval%*%c2
      tau_index = which(quadpts_eval == taunew)[1]
      x2_delay <- c(rep(x2_new[1], tau_index), D0quadbasismat_delay%*%c2)
      
      g1new <- 1/(1+(x2_delay/p_0new)^DDE_power) - mu_mnew*x1_new 
      #g2new <- x1 - mu_pnew*x2
      
      
      logacc_prob <- temperature*(-lambda*sum(quadwts_eval*((dx1_new-g1new)^2)) + lambda*sum(quadwts*((dx1-g1Old)^2)))
      if(runif(1, 0, 1) < exp(logacc_prob)){
        rt <- c(mu_mnew, mu_pnew, p_0new, taunew)
      }else{
        rt <- c(mu_mnew, mu_pnew, p_0new, tauOld)
      }  
      
      return(rt)
    }
    
    
    c2_ll <- function(xx2, y2, sigma2){
      rt <- -sum((y2 - xx2)^2/sigma2^2)
      return(rt)
    }
    
    
    c2_logprior <- function(lambda, parameters, dx1, dx2, x1, x2, x2_delay, quadwts, y2){
      mu_m <- parameters[1]
      mu_p <- parameters[2]
      p_0 <- parameters[3]
      tau <- parameters[4]
      
      g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
      g2 <- x1 - mu_p*x2
      
      rt <- -lambda*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))
      return(rt)
    }
    
    # Metropolis Hasting to sample basis coefficients c2
    MH_c2 <- function(c2_priorOld, c2_llOld, lambda, sigma2, parameters, dx1, x1, c2, quadwts, tau_index, basismat, D0quadbasismat, D1quadbasismat, D0quadbasismat_delay, y2, temperature, sigmac2){
      MHsigma <- 1
      newc2 <- c(rnorm(length(c2), c2, sigmac2))
      
      xx2 <- basismat%*%newc2
      dx2 <- D1quadbasismat%*%newc2
      x2 <- D0quadbasismat%*%newc2
      x2_delay <- c(rep(x2[1], tau_index), D0quadbasismat_delay%*%newc2)
      
      c2_priornew <- c2_logprior(lambda, parameters, dx1, dx2, x1, x2, x2_delay, quadwts, y2)
      c2_llnew <- c2_ll(xx2, y2, sigma2)
      logacc_prob <- temperature*(c2_priornew - c2_priorOld) + temperature*(c2_llnew - c2_llOld)
      if(runif(1, 0, 1) < exp(logacc_prob)){
        rt <- c(newc2)
      }else{
        rt <- c(c2)
      }
      return(list(c2_prior = c2_priornew, c2_ll = c2_llnew, c2 =rt)) 
    }
    
    # Gibbs step for measurement error sigma 
    GB_sigma <- function(y1, y2, x1, x2, temperature){
      temp1 <- rgamma(1, shape = 1+ length(y1)/2, rate = 1 + temperature*sum((y1 - x1)^2)/2)
      temp2 <- rgamma(1, shape = 1+ length(y2)/2, rate = 1 + temperature*sum((y2 - x2)^2)/2)
      return(c(sqrt(1/temp1), sqrt(1/temp2)))
    }
    
    # Gibbs move for lambda 
    GB_lambda <- function(alambda, blambda, parameters, M, dx1, dx2, x1, x2,  x2_delay, quadwts, temperature){
      mu_m <- parameters[1]
      mu_p <- parameters[2]
      p_0 <- parameters[3]
      tau <- parameters[4]
      
      g1 <- 1/(1+(x2_delay/p_0)^DDE_power) - mu_m*x1 
      g2 <- x1 - mu_p*x2
      
      b_lambda <- 1/(1/blambda+temperature*sum(quadwts*((dx1-g1)^2+(dx2-g2)^2))/2)
      return(rgamma(1, shape = alambda+temperature*M*2/2, scale = b_lambda))
    }
    
    
    # Gibbs move for c1 
    GB_c1 <- function(lambda, sigma1, parameters, dx1, dx2, x1, x2,  x2_delay, quadwts, D0quadbasismat, D1quadbasismat, y1, temperature, basismat){
      mu_m <- parameters[1]
      mu_p <- parameters[2]
      p_0 <- parameters[3]
      tau <- parameters[4]
      
      VV <- diag(quadwts)
      
      SigInv <- temperature*(lambda/2*(t(D0quadbasismat)%*%VV%*%D0quadbasismat+t(D1quadbasismat+mu_m*D0quadbasismat)%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+ t(basismat)%*%basismat/(2*sigma1^2))
      Sig <- solve(SigInv)
      mu <- t(temperature*(lambda/2*(t(dx2+mu_p*x2)%*%VV%*%D0quadbasismat+t(1/(1+(x2_delay/p_0)^DDE_power))%*%VV%*%(D1quadbasismat+mu_m*D0quadbasismat))+t(y1)%*%basismat/(2*sigma1^2))%*%Sig)
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
    
    
    output <- DE3(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  return(list(c = output$c, parameters = output$parameters, lambda = output$lambda, sigma = output$sigma, W = output$W))
}





