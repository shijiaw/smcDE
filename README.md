# smcDE

Summary
-------
smcDE is an R package of a recently proposed method  for estimating parameters in differential equations. The differential equations are solved via collocation, and inference are made by developed sequential Monte Carlo algorithm, based on constructed fully Bayesian scheme. 



Installation
------------
install_github(``shijiaw/smcDE'')

Input arguments
-----
- data: a list with each element being a trajectory with measurement error.
- times: a list with each element being the time we record data.
- seed: the random seeds for propagation and resampling in SMC.
- knots: the location we put knots.
- NP: the number of particles in SMC.
- CESSthresholds: relative conditional effective sample size (0, 1) in SMC.
- resampleThreshold: the threshold triggering resampling in SMC.
- alambda: hyper-parameter for tuning parameter lambda.
- blambda: hyper-parameter for tuning parameter lambda.
- sigmac: the standard deviation of the basis coefficients c.
- DEmodel: options for DE models: 1 (ODE), 2 (DDE1), 3 (DDE2) in manuscript.