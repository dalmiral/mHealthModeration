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

attributes(rsnmm(5, 5, beta0 = c(-0.8, 0, 0.2, 0)))

## low moderation
omit1 <- sim(n, tmax, M, scenario = "omit", beta0 = c(-0.8, 0, 0.2, 0))

## medium moderation
omit2 <- sim(n, tmax, M, scenario = "omit", beta0 = c(-0.8, 0, 0.5, 0))

## high moderation
omit3 <- sim(n, tmax, M, scenario = "omit", beta0 = c(-0.8, 0, 0.8, 0))

save(omit1, omit2, omit3, file = "sim-omit.RData")
