library("foreach")
library("doParallel")
library("parallel")
source("init.R")
source("sim.R")

## set number of Monte Carlo replicates
M <- 1000

## set number of threads to use for parallel processing and the random seed
## (nb: these two values ensure that the results are replicable)
cores <- 4
seed <- 0

cl <- makeCluster(getOption("cl.cores", cores))
clusterEvalQ(cl, source("init.R"))
registerDoParallel(cl)

sim.omit <- function() {
  out <- NULL
  ## low, medium and high degrees of moderation by state
  for (b in c(0.2, 0.5, 0.8)) {
    for (n in c(30, 60)) {
      for (tmax in c(30, 50)) {
        clusterSetRNGStream(cl, seed)
        out <-
          rbind(out,
                cbind(level = paste("$\\beta_{11}^* = ", b, "$", sep = ""),
                      sim(n, tmax, M,
                          ## regress response on state and proximal treatment,
                          ## ignoring the underlying interaction between the two
                          y.formula = list(w = y ~ state + I(a - pn),
                                           u.ind = y ~ state + a,
                                           u.ar1 = y ~ state + a,
                                           u.exch = y ~ state + a),
                          y.names = c(w = "Weighted and centered",
                                      u.ind = "GEE independence",
                                      u.ar1 = "GEE AR(1)",
                                      u.exch = "GEE exchangeable"),
                          ## term labels for proximal treatment
                          y.label = list(w = "I(a - pn)",
                                         u.ind = "a", u.ar1 = "a", u.exch = "a"),
                          ## specify weights and working correlation structure
                          y.args = list(w = list(wn = "pn", wd = "prob"),
                                        u.ind = list(),
                                        u.ar1 = list(corstr = "ar1"),
                                        u.exch = list(corstr = "exch")),
                          ## specify weight numerator model
                          a.formula = list(pn = a ~ 1),
                          a.names = c(pn = "Intercept-only"),
                          ## use default generative model, but with the specified
                          ## level of moderation by the time-varying state
                          beta0 = c(-0.2, 0, 0, b, 0))))
      }
    }
  }
  out
}

omit <- sim.omit()
save(omit, file = "sim-omit.RData")

stopCluster(cl)
