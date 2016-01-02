source("utils.R")
library("foreach")

sim <- function(n = 25, tmax = 25, M = 1000, cp = 0.95, proximal = TRUE,
                scenario = c("omit", "stable", "ar1", "truth"),
                effects = TRUE, control = rsnmm.control(), ...)
{
  if (missing(control))
    control <- rsnmm.control(...)
  else
    control <- do.call("rsnmm.control", control)
  scenario <- match.arg(scenario)
  d <- rsnmm(n = n, tmax = tmax, control = control)
  at <- attributes(d)
  tval <- (at$lag + !proximal):max(d$time)
  lag.prefix <- c("lag", "")[proximal + 1]
  ## omit these terms from treatment (main) effect, regardless of the truth
  o <- paste0(lag.prefix, "vary")
  truth <- if (proximal) params(at$beta0, omit = o)
           else params(at$beta1, omit = o)
  ## lists of quoted variables names for cbind, via do.call
  ## - weight "stabilizers"
  ht <- qnames(params(at$eta, omit = "laga"))
  ## - treatment effect
  ha <- qnames(truth, prefix = paste0(lag.prefix, "a"), omit = o)
  hac <- qnames(truth, prefix = paste0(lag.prefix, "a.center"), omit = o)
  ## - main effect
  if (scenario == "omit")
    o <- ""
  else if (scenario == "ar1")
    k <- c("one", "lagy")
  else
    k <- ""
  hm <- if (proximal) with(at, params(c(mu, xi0), keep = k))
        else with(at, params(c(mu, xi1)))
  hm <- hmc <- unique(c(qnames(c(truth, hm), keep = k, omit = o)))
  if (proximal & any(at$beta1 != 0)) {
    hm <- unique(c(hm, qnames(params(at$beta1, omit = o), prefix = "laga")))
    hmc <- unique(c(hmc, qnames(params(at$beta1, omit = o),
                                prefix = "laga.center")))
  }
  cat("\nGenerative model parameters\n\n")
  print(data.frame(t(at$beta0), row.names = "  Proximal effect"))
  cat("\n")
  print(data.frame(t(at$beta1), row.names = "  Delayed effect"))
  cat("\n")
  print(data.frame(t(with(at, c(mu, xi0, xi1))), row.names = "  Main effect"))
  cat("\n")
  print(data.frame(t(at$coef.vary), row.names = "  Time-varying covariate"))
  cat("\n")
  print(data.frame(t(at$eta), row.names = "  Treatment"))
  cat("\nAnalysis model terms - weighting or routine")
  cat("\n\n  Treatment:", paste(ha, collapse = ", "))
  cat("\n  Main:", paste0(hm, collapse = ", "))
  if (scenario != "stable") {
    cat("\n\nAnalysis model terms - centering")
    cat("\n\n  Treatment:", paste0(hac, collapse = ", "))
    cat("\n  Main:", paste0(hmc, collapse = ", "))
  }
  cat("\n\n")
  ## model fitter via gack
  ## d is the data frame for the replicate
  fit.gee <- function(mod = ha, main = hm, wgt = quote(w), gn = NULL,
                      corstr = c("independence", "fixed"), name = "Weighted",
                      ...) {
    corstr <- match.arg(corstr)
    k <- length(mod)
    if (corstr == "independence") {
      corr <- diag(length(tval))
      zcor <- NULL
    }
    else {
      corr <- attributes(d)$corr[1:length(tval), 1:length(tval)]
      zcor <- fixed2Zcor(corr, rep(1:n, each = length(tval)),
                         rep(1:length(tval), n))
    }
    f <- try(with(subset(d, time %in% tval),
                  geese.glm(do.call("cbind", c(mod, main)), y = y, id = id,
                            weights = eval(wgt), corstr = corstr, zcor = zcor,
                            ...)))
    if (inherits(f, "try-error"))
      stop(name, " response model failed at (n, T, i) = (", n, ", ", tmax, ", ",
           i, ")\n", sep = "")
    f$trtcoef <- rep(0, length(f$geese$beta))
    f$trtcoef[1:k] <- 1
    v <- working.covariance(f, invert = TRUE, wcor = corr)
    b <- bread.geeglm(f, wcovinv = v)
    m <- meat.geeglm(f, gn = gn, wcovinv = v, small = n < 50, lag = !proximal)
    list(name = name, k = k, beta = f$geese$beta[1:k],
         se = sqrt(diag(f$geese$vbeta))[1:k],
         sec = sqrt(diag(b %*% m %*% t(b)))[1:k])
  }
  ## for each Monte Carlo replicate
  res <- foreach(i = 1:M, .combine = rbind) %dopar% {
    d <- rsnmm(n = n, tmax = tmax, control = control)
    ## calculate weights
    fitn <- if (scenario != "ar1")
              try(with(subset(d, time %in% tval),
                       geese.glm(x = cbind(one), y = a, id = id,
                                 family = binomial())), silent = TRUE)
            else NULL
    if (inherits(fitn, "try-error"))
      stop("Overall probability estimate failed at (n, T, i) = (", n, ", ",
           tmax, ", ", i, ")\n", sep = "")
    d$rho <- if (scenario != "ar1") fitn$fitted.values[1]
             else 0.5
    d$w <- with(d, (rho / prob)^a * ((1 - rho) / (1 - prob))^(1 - a))
    if (!proximal)
      d$w <- with(d, delay(id, time, w))
    if (scenario == "stable") {
      fitn.vary <- try(with(subset(d, time %in% tval),
                            geese.glm(x = do.call("cbind", ht), y = a,
                                      id = id, family = binomial())),
                       silent = TRUE)
      if (inherits(fitn.vary, "try-error"))
        stop("Time-varying probability estimate failed at (n, T, i) = (",
             n, ", ", tmax, ", ", i, ")\n", sep = "")
      d$rho.vary <- 0
      d$rho.vary[d$time %in% tval] <- fitn.vary$fitted.values
      d$w.vary <- with(d, (rho.vary / prob)^a
                       * ((1 - rho.vary) / (1 - prob))^(1 - a))
      if (!proximal)
        d$w.vary <- with(d, delay(id, time, w.vary))
    }
    ## fit response models
    fitc <- fitr <- fitwa <- fitra <- fitca <- NULL
    fitw <- fit.gee(gn = fitn)
    if (scenario == "stable")
      fitwa <- fit.gee(wgt = quote(w.vary), gn = fitn.vary,
                       name = "Weighted time-varying stabilizer")
    else if (scenario == "ar1")
      fitwa <- fit.gee(gn = fitn, corstr = "fixed",
                       name = "Weighted fixed AR(1)")
    else {
      fitr <- fit.gee(wgt = quote(one), name = "Routine")
      fitc <- fit.gee(hac, hmc, wgt = quote(one), name = "Centered")
    }
    fits <- list(w = fitw, wa = fitwa, c = fitc, r = fitr)
    fits <- fits[!sapply(fits, is.null)]
    ds <- data.frame(truth, names(truth), 1, min(with(d, prob[time %in% tval])),
                     max(with(d, prob[time %in% tval])))
    ds <-
      cbind(ds,
            do.call("cbind", lapply(fits, function(f) with(f, rep(name, k)))),
            do.call("cbind", lapply(fits, function(f) with(f, beta))),
            do.call("cbind", lapply(fits, function(f) with(f, se))),
            do.call("cbind", lapply(fits, function(f) with(f, sec))),
            do.call("cbind", lapply(fits, function(f)
             with(f, inside.ci(ds$truth, beta, se, cp)))),
            do.call("cbind", lapply(fits, function(f)
              with(f, inside.ci(ds$truth, beta, sec, cp, qt, df = n - k)))))
    ds <- setNames(ds, c("truth", "term", "nrep", "pmin", "pmax",
                         paste("fit", names(fits), sep = "."),
                         paste("est", names(fits), sep = "."),
                         paste("se", names(fits), sep = "."),
                         paste("sec", names(fits), sep = "."),
                         paste("cp", names(fits), sep = "."),
                         paste("cpc", names(fits), sep = ".")))
    ds
  }
  ifit <- c(1:2, grep("^fit\\.", names(res)))
  iavg <- c(1:2, grep("^(est|se|sec|cp|cpc)\\.", names(res)))
  isd <- c(1:2, grep("^est\\.", names(res)))
  data.frame(n = n, tmax = tmax, lag = 1 + !proximal,
             res[1, ifit, drop = FALSE],
             aggregate(nrep ~ term + truth,
                       data = res, FUN = sum)[, 3, drop = FALSE],
             aggregate(pmin ~ term + truth,
                       data = res, FUN = min)[, 3, drop = FALSE],
             aggregate(pmax ~ term + truth,
                       data = res, FUN = max)[, 3, drop = FALSE],
             aggregate(. ~ term + truth,
                       data = res[, iavg], FUN = mean)[, -(1:2)],
             setNames(aggregate(. ~ term + truth,
                                data = res[, isd], FUN = sd)[, -(1:2)],
                      gsub("^est\\.", "sd.", names(res[, isd[-(1:2)]]))),
             row.names = NULL)
}
