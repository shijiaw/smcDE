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
    source("DE1.R")
    output <- DE1(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  if(DEmodel == 2){
    source("DE2.R")
    output <- DE2(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  if(DEmodel == 3){
    source(("DE3.R"))
    output <- DE3(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac)
  }
  
  return(list(c = output$c, parameters = output$parameters, lambda = output$lambda, sigma = output$sigma, W = output$W))
}





