library("foreach")
library("doParallel")
library("parallel")
source("init.R")
source("sim.R")

## set number of Monte Carlo replicates
M <- 1000

## set number of threads to use for parallel processing and the random seed
## (nb: these two values ensure that the results are replicable)
cores <- 2
seed <- 0

cl <- makeCluster(getOption("cl.cores", cores))
clusterEvalQ(cl, source("init.R"))
registerDoParallel(cl)

## control parameters for the generative model
beta0 <- c(-0.2, 0, 0, 0, 0)     # unmoderated proximal effect
beta1 <- c(-0.1, 0, 0, 0)        # unmoderated delayed effect
coef.state <- c(0, 0, 0, 0, 0.1) # state depends on past treatment
eta <- rep(0, 5)                 # treatment probability = 1/2

sim.ar1 <- function() {
  out <- NULL
  for (n in c(30, 60)) {
    for (tmax in c(30, 50)) {
      ## obtain true correlation matrix, trimmed down to effective size
      ## ("effective" observations avoid (lags of) initial values)
      attrib <- attributes(rsnmm(n, tmax, beta0 = beta0, beta1 = beta1,
                                 coef.state = coef.state, eta = eta))
      ## by default the true correlation structure is AR(1) with (u, t)th error
      ## correlation sqrt(0.5)^abs(u - t)
      cormatrix <- attrib$cormatrix[1:(tmax - attrib$lag + 1),
                                    1:(tmax - attrib$lag + 1)]
      clusterSetRNGStream(cl, seed)
      out <- rbind(out,
                   sim(n, tmax, M,
                       ## regress response on proximal treatment, centered by the
                       ## true treatment probability
                       y.formula = list(indep = y ~ state + I(a - prob),
                                        ar1 = y ~ state + I(a - prob)),
                       y.names = c(indep = "Independence",
                                   ar1 = "AR(1)"),
                       y.label = list(indep = "I(a - prob)",
                                      ar1 = "I(a - prob)"),
                       ## employ different working correlation structures
                       y.args = list(indep = list(),
                                     ar1 = list(corstr = "userdefined",
                                                wcor = cormatrix)),
                       a.formula = NULL, a.names = NULL,
                       beta0 = beta0, beta1 = beta1, coef.state = coef.state,
                       eta = eta))
    }
  }
  out
}

ar1 <- sim.ar1()
save(ar1, file = "sim-ar1.RData")

stopCluster(cl)
