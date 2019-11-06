# smcDE
An R package described in our manuscript ``Adaptive semiparametric Bayesian differential equations via sequential Monte Carlo methods''

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