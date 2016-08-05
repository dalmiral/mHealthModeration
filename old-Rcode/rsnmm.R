rsnmm <- function(n, tmax, control, ...)
{
  if (missing(control))
    control <- rsnmm.control(...)
  else
    control <- do.call("rsnmm.control", control)
  tmax <- tmax + (tmax %% 2) + 1
  v <- matrix(0, tmax, tmax)
  control$alpha <- control$corr
  control$corr <- matrix(control$corr^(abs(row(v) - col(v))), tmax, tmax)
  v <- with(control, corr * sd[1]^2)
  err <- apply(matrix(rnorm(n * tmax), tmax), 2, function(x) t(chol(v) %*% x))
  d <- .C("rsnmm",
          n = as.integer(n),
          tmax = as.integer(tmax),
          beta = with(control, as.double(c(beta0, beta1))),
          eta = as.double(control$eta),
          mu = as.double(control$mu),
          xi = with(control, as.double(c(xi0, xi1))),
          coef.vary = as.double(control$coef.vary),
          base = as.double(rep(rnorm(n), each = tmax)),
          vary = as.integer(rep(0, n * tmax)),
          a = as.integer(rep(0, n * tmax)),
          prob = as.double(rep(0, n * tmax)),
          y = as.double(rep(0, n * tmax)),
          err = as.double(err),
          vary.center = as.double(rep(0, n*tmax)),
          a.center = as.double(rep(0, n*tmax)),
          sd = as.double(control$sd))
  d <- data.frame(id = rep(1:n, each = tmax), time = rep(0:(tmax - 1), n),
                  base = d$base, vary = d$vary, a = d$a, y = d$y, err = d$err,
                  prob = d$p, a.center = d$a.center, vary.center = d$vary.center,
                  one = 1)
  d$lagy <- with(d, delay(id, time, y))
  d$lag2y <- with(d, delay(id, time, y, 2))
  d$lagerr <- with(d, delay(id, time, err))
  d$laga <- with(d, delay(id, time, a))
  d$lag2a <- with(d, delay(id, time, a, 2))
  d$lagprob <- with(d, delay(id, time, prob))
  d$lag2prob <- with(d, delay(id, time, prob, 2))
  d$laga.center <- with(d, delay(id, time, a.center))
  d$lag2a.center <- with(d, delay(id, time, a.center, 2))
  d$lagvary <- with(d, delay(id, time, vary))
  d$lagvary.center <- with(d, delay(id, time, vary.center))
  rownames(d) <- NULL
  attributes(d) <- c(attributes(d), control)
  d
}

rsnmm.control <- function(origin = 1, sd = 1, corr = 0.5,
                          beta0 = c(-0.8, 0, 0.8, 0), beta1 = rep(0, 3),
                          eta = c(0, 0, 0.8, -0.8, 0), mu = c(0, 0),
                          xi0 = c(0.8, 0), xi1 = 0, coef.vary = rep(0, 4),
                          lag = 3 + any(beta1 != 0))
{
  list(origin = 1, sd = sd, corr = corr, lag = lag,
       beta0 = setNames(beta0, c("one", "base", "vary", "laga")),
       beta1 = setNames(beta1, c("one", "base", "lagvary")),
       eta = setNames(eta, c("one", "base", "vary", "laga", "lagy")),
       mu = setNames(mu, c("one", "base")),
       xi0 = setNames(xi0, c("vary", "lagy")), xi1 = setNames(xi1, "lagvary"),
       coef.vary = setNames(coef.vary, c("one", "base", "lagvary", "laga")))
}
