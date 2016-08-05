library("foreach")
library("doParallel")
library("parallel")
source("sim.R")

cl <- makeCluster(getOption("cl.cores", 4))
clusterSetRNGStream(cl, 0)
clusterEvalQ(cl, source("utils.R"))
registerDoParallel(cl)

n <- 100
tmax <- 100
M <- 1000

attributes(rsnmm(5, 5, beta0 = c(-0.8, 0, 0, 0)))

stable <- sim(n, tmax, M, scenario = "stable", beta0 = c(-0.8, 0, 0, 0))

save(stable, file = "sim-stable.RData")
