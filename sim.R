sim <- function(n = 30, tmax = 30, M = 1000,
                ## response regression models
                y.formula = list(w = y ~ I(a - pn) * (base + state),
                                 u = y ~ a * (base + state)),
                ## names for each regression model
                y.names = c(w = "Weighted and centered",
                            u = "GEE AR(1)"),
                ## labels for regression terms of the treatment effect
                y.label = list(w = "I(a - pn)"),
                ## names of the treatment probability models or variables used
                ## for the weight numerator ('wn') or denominator ('wd') and
                ## arguments for the estimation routine
                y.args = list(w = list(wn = "pn", wd = "pd"),
                              u = list(corstr = "ar1")),
                ## treatment probability models named in 'y.args'
                a.formula = list(pn = a ~ lag1a,
                                 pd = a ~ lag1a + state),
                ## names for each treatment probability model
                a.names = c(pn = "Last treatment",
                            pd = "Last treatment and current state"),
                ## proximal (0) or delayed (1) treatment effect?
                lag = 0,
                ## print generative and analysis model details
                verbose = TRUE,
                ## control parameters for 'rsnmm'
                control, ...) {
  control <- if (missing(control)) rsnmm.control(...)
             else control <- do.call("rsnmm.control", control)
  ## times to use in the model fit
  runin.fita <- control$lag
  runin.fity <- control$lag + lag
  ## retrieve causal control parameter values
  ## nb: if the regression models 'y.formula' average over an underlying
  ##     moderator these will not represent the true causal effect unless this
  ##     moderator has conditional mean zero
  y.coef <- mapply(which.terms, x = y.formula, label = y.label,
                   stripnames = TRUE, SIMPLIFY = FALSE)
  truth <- control[[paste0("beta", lag)]]
  truth <- truth[Reduce("intersect", lapply(y.coef, names))]
  y.coef <- lapply(y.coef, function(x) x[names(truth)])
  ## corresponding treatment probability models
  ## nb: we avoid delayed evaluation in 'y.args' (e.g. passing a 'weights'
  ##     argument directly) to avoid scoping issues in 'foreach'
  if (!is.null(a.formula)) {
    y.prob <- lapply(y.args, function(x) do.call("c", x[c("wn", "wd")]))
    y.prob <- lapply(y.prob, function(x) x[x %in% names(a.formula)])
  }
  else y.prob <- lapply(y.formula, function(x) list())
  ## print generative and analysis model properties
  if (verbose) {
    cat("\nGenerative model attributes\n\n")
    print(control)
    cat("Analysis models\n\n")
    mapply(function(f, nm) write.table(cbind("  ", nm, ": y ~ ",
                                            as.character(f)[3]), sep = "",
                                      row.names = FALSE, col.names = FALSE,
                                      quote = FALSE, eol = "\n\n"),
           f = y.formula, nm = y.names)
    cat("Treatment probability models\n\n")
    mapply(function(f, nm) write.table(cbind("  ", nm, ": a ~ ",
                                            as.character(f)[3]), sep = "",
                                      row.names = FALSE, col.names = FALSE,
                                      quote = FALSE, eol = "\n\n"),
           f = a.formula, nm = a.names)
  }
  ## general model fitter
  ## nb: d is the data frame for the replicate
  fitter <- function(formula, args, prob, coef, label, response = "y",
                     addvar = NULL) {
    if (response == "a") {
      args$family <- binomial()
      runin <- runin.fita
    }
    else runin <- runin.fity
    r <- which(d$time >= runin)
    l <- list(x = model.matrix(formula, data = d[r, ]), y = d[r, response])
    if (is.null(args$wn) & is.null(args$wn)) l$w <- rep(1, nrow(d))
    else {
      l$w <- ifelse(d[, "a"] == 1, d[, args$wn] / d[, args$wd],
                    (1 - d[, args$wn]) / (1 - d[, args$wd]))
      args[c("wn", "wd")] <- NULL
    }
    l$w <- l$w * d$avail
    if (lag) l$w <- delay(d$id, d$time, l$w, lag)
    l$w <- l$w[r]
    if (!is.null(args$corstr)) {
      fun <- "geese.glm"
      l$id <- d$id[r]
    }
    else if (!is.null(args$family)) fun <- "glm.fit"
    else fun <- "lm.wfit"
    fit <- do.call(fun, c(l, args))
    if (!inherits(fit, "geeglm")) {
      fit <- glm2gee(fit, d$id[r])
      fit$terms <- terms(formula)
      fit$geese$X <- l$x
      fit$y <- l$y
    }
    if (!is.null(addvar)) {
      newvar <- paste0(c("", "lag1"), addvar)
      d[, newvar] <<- NA
      d[r, newvar[1]] <<- fit$fitted.values
      d[, newvar[2]] <<- delay(d$id, d$time, d[, newvar[1]])
    }
    else {
      ## usual variance sandwich estimator
      fit$vcov <- vcov.geeglm(fit)
      est <- estimate(fit)[coef, 1:4, drop = FALSE]
      ## correction for any estimates in weights
      l <- if (length(prob)) setNames(fita[prob], gsub("^w", "p", names(prob)))
           else NULL
      fit$vcov <- NULL
      fit$vcov <- do.call("vcov.geeglm", c(list(x = fit, label = label), l))
      estc <- estimate(fit)[coef, 1:4, drop = FALSE]
      fit <- data.frame(moderator = names(coef), truth = truth,
                        est = est[, "Estimate"], se = est[, "SE"],
                        lcl = est[, "95% LCL"], ucl = est[, "95% UCL"],
                        sec = estc[, "SE"], lclc = estc[, "95% LCL"],
                        uclc = estc[, "95% UCL"], row.names = NULL)
    }
    fit
  }
  fita <- list()
  ## for each replicate...
  out <- foreach(m = 1:M, .combine = "rbind") %dopar% {
    ## ... generate data
    d <- rsnmm(n, tmax, control = control)
    d$pn <- d$pd <- d$prob
    ## ... fit treatment probability models
    if (!is.null(a.formula))
      fita <- mapply(fitter, formula = a.formula, addvar = names(a.formula),
                     MoreArgs = list(args = list(), prob = list(),
                                     coef = list(), label = list(),
                                     response = "a"), SIMPLIFY = FALSE)
    ## ... fit response models
    fity <- mapply(fitter, formula = y.formula, args = y.args, prob = y.prob,
                   coef = y.coef, label = y.label, SIMPLIFY = FALSE)
    fity <- mapply(function(nm, d) data.frame(iter = m, method = nm, d,
                                              row.names = NULL),
                   nm = y.names[names(fity)], d = fity, SIMPLIFY = FALSE)
    out <- do.call("rbind", setNames(fity, NULL))
  }
  out <- data.frame(n, tmax, out)
  ## 95% CI coverage probability using uncorrected SEs
  out$cp <- with(out, lcl <= truth & truth <= ucl)
  ## coverage probability using SEs corrected for estimates in weights
  out$cpc <- with(out, lclc <= truth & truth <= uclc)
  ## root MSE
  out$rmse <- with(out, (est - truth)^2)
  ## mean and SD estimate, number of replicates
  out <- cbind(aggregate(cbind(est, se, sec, cp, cpc, rmse) ~
                           method + moderator + truth + n + tmax,
                         data = out, FUN = mean),
               sd = aggregate(est ~ method + moderator + truth + n + tmax,
                              data = out, FUN = sd)$est,
               iter = aggregate(iter ~ method + moderator + truth + n + tmax,
                                data = out,
                                FUN = function(x) length(unique(x)))$iter)
  out$rmse <- sqrt(out$rmse)
  out
}
