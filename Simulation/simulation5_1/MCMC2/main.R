rm(list=ls())

nSamples <- 400000
filename_theta1 = './theta1'
filename_theta2 = './theta2'
filename_sigma1 = './sigma1'
filename_sigma2 = './sigma2'
filename_initial1 = './initial1'
filename_initial2 = './initial2'
source("simulator.R")
source("MetropolisHasting.R")

thin_iter <- (0:9999)*40+1

theta1_thin = thetaStore[thin_iter, 1]
theta2_thin = thetaStore[thin_iter, 2]
sigma1_thin = sigmaStore[thin_iter, 1]
sigma2_thin = sigmaStore[thin_iter, 2]
initial1_thin = initialStore[thin_iter, 1]
initial2_thin = initialStore[thin_iter, 2]

write.table(matrix(round(theta1_thin, digits = 5),nr=1), file = paste(filename_theta1,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(theta2_thin, digits = 5),nr=1), file = paste(filename_theta2,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(sigma1_thin, digits = 5),nr=1), file = paste(filename_sigma1,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(sigma2_thin, digits = 5),nr=1), file = paste(filename_sigma2,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(initial1_thin, digits = 5),nr=1), file = paste(filename_initial1,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)
write.table(matrix(round(initial2_thin, digits = 5),nr=1), file = paste(filename_initial2,".txt",sep=""),  row.names = FALSE,
            col.names = FALSE)


