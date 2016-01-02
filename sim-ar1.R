library("foreach")
library("doParallel")
library("parallel")
source("sim.R")

cl <- makeCluster(getOption("cl.cores", 4))
clusterSetRNGStream(cl, 0)
clusterEvalQ(cl, source("utils.R"))
registerDoParallel(cl)

n <- 50
tmax <- 100
M <- 1000

attributes(rsnmm(5, 5, beta0 = c(-0.8, 0, 0, 0), eta = rep(0, 5),
                 xi0 = c(0.8, 0), coef.vary = c(0, 0, 0, 0.8)))

ar1 <- sim(n, tmax, M, scenario = "ar1", beta0 = c(-0.8, 0, 0, 0),
           eta = rep(0, 5), xi0 = c(0.8, 0), coef.vary = c(0, 0, 0, 0.8))

save(ar1, file = "sim-ar1.RData")
