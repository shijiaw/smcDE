# smcDE

Summary
-------
smcDE is an R package of a recently proposed method  for estimating parameters in differential equations. The differential equations are solved via collocation, and inference are made by developed sequential Monte Carlo algorithm, based on constructed fully Bayesian scheme. 



Installation
------------
install_github(``shijiaw/smcDE'')

Quick Start
------------
smcDE(data, times, seed, knots, CESSthresholds, NP, resampleThreshold, alambda, blambda, sigmac, DEmodel)

Input Arguments
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

Output
-----
- c: a list of particles for basis coefficient (c1 and c2).
- parameter: particles for model parameters.
- lambda: particles for tuning parameter lambda. 
- sigma: particles for sigma.
- W: normalized weights for particles.

Demo
-----
We refer users to Folder `Demos' for  examples. In `Demos', folder `Demo1' includes an example for ODE parameter estimation; folder `Demo2' includes an example for our first DDE example; folder `Demo3' includes an example for our second DDE example. 

Simulation
-----
The folder `Simulation' includes setups for simulation studies in manuscript. The results can be reproduced by running `main.R'.

