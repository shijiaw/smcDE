rm(list=ls())

install.packages("devtools")
library("devtools")

# install package
install_github("shijiaw/smcDE")
library("smcDE")

# load data
load("data_DDE2.Rdata")
load("time.Rdata")

# run SMC
output = smcDE(data = data_DDE2, times = time, seed = 364, knots = seq(time[[1]][1], max(time[[1]]),by = 20), CESSthresholds = 0.8, NP = 300, resampleThreshold = 0.5, alambda = 1, blambda = 1, sigmac = 100, DEmodel = 3)



