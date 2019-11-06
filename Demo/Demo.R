rm(list=ls())

install.packages("devtools")
library("devtools")

# install package
install_github("shijiaw/smcDE")
library("smcDE")

# load data
load("data1.Rdata")
load("time1.Rdata")

# run SMC
output = smcDE(data = data, times = time, seed = 1, knots = seq(time[[1]][1], max(time[[1]]),by = 5), CESSthresholds = 0.8, NP = 500, resampleThreshold = 0.5, alambda = 1, blambda = 1, sigmac = 100, DEmodel = 1)



