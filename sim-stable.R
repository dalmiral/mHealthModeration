library("foreach")
library("doParallel")
library("parallel")
source("init.R")
source("sim.R")

## set number of Monte Carlo replicates
M <- 1000

## set number of threads to use for parallel processing and the random seed
## nb: these two values ensure that the results are replicable
cores <- 2
seed <- 0

cl <- makeCluster(getOption("cl.cores", cores))
clusterEvalQ(cl, source("init.R"))
registerDoParallel(cl)

sim.stable <- function() {
  out <- NULL
  for (n in c(30,60)) {
    for (tmax in c(30,50)) {
      clusterSetRNGStream(cl, seed)
      out <- rbind(out,
                   sim(n, tmax, M,
                       ## regress response on proximal treatment, centered by a
                       ## probability that is either (i) constant over time or
                       ## (ii) time-varying
                       y.formula = list(fixed = y ~ state + I(a - pfixed),
                                        vary = y ~ state + I(a - pvary)),
                       y.names = c(fixed = "Constant in $t$ (i)",
                                   vary = "Depends on $S_t$ (ii)"),
                       y.label = list(fixed = "I(a - pfixed)",
                                      vary = "I(a - pvary)"),
                       ## weight regression using the true treatment probability
                       ## in the denominator
                       y.args = list(fixed = list(wn = "pfixed", wd = "prob"),
                                     vary = list(wn = "pvary", wd = "prob")),
                       ## model numerator probability with (i) intercept only and
                       ## (ii) state, which is time-varying
                       a.formula = list(pfixed = a ~ 1,
                                        pvary = a ~ state),
                       a.names = c(pfixed = "Constant in $t$ (i)",
                                   pvary = "Depends on $S_t$ (ii)"),
                       ## use default generative model, but with medium level of
                       ## moderation by state and an unmoderated delayed effect
                       beta0 = c(-0.2, 0, 0, 0.5, 0),
                       beta1 = c(-0.1, 0, 0, 0)
                       ))
    }
  }
  out
}

stable <- sim.stable()
save(stable, file = "sim-stable.RData")

stopCluster(cl)
