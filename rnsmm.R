rsnmm <- function(n, tmax, control, ...) {
  control <- if (missing(control)) rsnmm.control(...)
             else do.call("rsnmm.control", control)
  tmax <- tmax + (tmax %% 2) + 1
  time <- rep(0:(tmax - 1), n)
  tfun <- do.call("data.frame", lapply(control$tfun, function(f) f(time, tmax)))
  ## provisional error
  err <- rnorm(n * tmax, sd = control$sd)
  coef.err <- 0
  control$cormatrix <- matrix(control$coralpha, tmax, tmax)
  if (control$corstr == "exchangeable") {
    diag(control$cormatrix) <- 1
    ## add cluster error component
    err <- control$coralpha * rep(rnorm(n, sd = 1 / control$coralpha),
                                  each = tmax) + err
  }
  else if (control$corstr == "ar1") {
    coef.err <- control$coralpha
    err[time > 0] <- rnorm(n * (tmax - 1), sd = sqrt(1 - control$coralpha^2))
    control$cormatrix <- matrix(with(control,
                                     coralpha^(abs(row(cormatrix) -
                                                   col(cormatrix)))), tmax, tmax)
  }
  expit <- function(a) exp(a) / (1 + exp(a))
  avail <- avail.center <- base <-
    state <- state.center <- a <- a.center <- prob <- y <- rep(0, n * tmax)
  for (u in 2:tmax) {
    k <- which(time == (u - 1))
    r <- expit(control$coef.avail[1] +
               control$coef.avail[2] * tfun$tavail[u] +
               control$coef.avail[3] * a[k - 1] +
               control$coef.avail[4] * y[k - 1])
    avail[k] <- sapply(r, rbinom, n = 1, size = 1)
    avail.center[k] <- avail[k] - r
    q <- expit(control$coef.state[1] +
               control$coef.state[2] * tfun$tstate[u] +
               control$coef.state[3] * base[k - 1] +
               control$coef.state[4] * state[k - 1] +
               control$coef.state[4] * a[k - 1])
    state[k] <- sapply(q, rbinom, n = 1, size = 1)
    state.center[k] <- state[k] - q
    prob[k] <- avail[k] * expit(control$eta[1] +
                                control$eta[2] * base[k] +
                                control$eta[3] * state[k] +
                                control$eta[4] * a[k - 1] +
                                control$eta[5] * y[k - 1])
    a[k] <- sapply(prob[k], rbinom, n = 1, size = 1)
    a.center[k] <- a[k] - prob[k]
    yc <-
      control$mu[1] +
      control$mu[2] * tfun$ty[u] +
      control$mu[3] * base[k] +
      a.center[k] * (control$beta0[1] +
                     control$beta0[2] * tfun$tmod[u] +
                     control$beta0[3] * base[k] +
                     control$beta0[4] * state[k] +
                     control$beta0[5] * a[k - 1]) +
      a.center[k - 1] * (control$beta1[1] +
                         control$beta1[2] * tfun$tmod[u - 1] +
                         control$beta1[3] * base[k - 1] +
                         control$beta1[4] * state[k - 1]) +
      control$theta0[1] * avail.center[k] +
      control$theta0[2] * state.center[k] +
      control$theta1[1] * avail.center[k - 1] +
      control$theta1[2] * state.center[k - 1]
    err[k] <- err[k] + coef.err * err[k - 1]
    y[k] <- yc + err[k]
  }
  d <- data.frame(id = rep(1:n, each = tmax), time,
                  ty = tfun$ty, tmod = tfun$tmod, tavail = tfun$tavail,
                  tstate = tfun$tstate, base, state, a, y, err,
                  avail, prob, a.center, state.center, avail.center)
  ## nb: for a given row, y is the proximal responsel
  d$lag1y <- with(d, delay(id, time, y))
  d$lag2y <- with(d, delay(id, time, y, 2))
  d$lag1err <- with(d, delay(id, time, err))
  d$lag1avail <- with(d, delay(id, time, avail))
  d$lag1avail.center <- with(d, delay(id, time, avail.center))
  d$lag2avail <- with(d, delay(id, time, avail, 2))
  d$lag2avail.center <- with(d, delay(id, time, avail.center, 2))
  d$lag1a <- with(d, delay(id, time, a))
  d$lag2a <- with(d, delay(id, time, a, 2))
  d$lag1prob <- with(d, delay(id, time, prob))
  d$lag2prob <- with(d, delay(id, time, prob, 2))
  d$lag1a.center <- with(d, delay(id, time, a.center))
  d$lag2a.center <- with(d, delay(id, time, a.center, 2))
  d$lag1tmod <- with(d, delay(id, time, tmod))
  d$lag2tmod <- with(d, delay(id, time, tmod, 2))
  d$lag1state <- with(d, delay(id, time, state))
  d$lag1state.center <- with(d, delay(id, time, state.center))
  rownames(d) <- NULL
  attributes(d) <- c(attributes(d), control)
  d
}

rsnmm.control <- function(origin = 1, sd = 1,
                          coralpha = 0.5,
                          corstr = c("independence", "ar1", "exchangeable"),
                          beta0 = c(-0.2, 0, 0, 0.2, 0), beta1 = rep(0, 4),
                          eta = c(0, 0, 0.8, -0.8, 0), mu = rep(0, 3),
                          theta0 = c(0, 0.8), theta1 = c(0, 0),
                          coef.avail = c(100, rep(0, 3)), coef.state = rep(0, 5),
                          tfun = NULL, lag = 3 + any(beta1 != 0)) {
  corstr <- match.arg(corstr)
  if (is.null(tfun))
    tfun <- rep(list(function(tcur, tmax) rep(0, length(tcur))), 4)
  list(origin = 1, lag = lag,
       ## error SD, correlation - AR(1) simplified by assuming unit variance
       sd = sd^(corstr != "ar1"),
       coralpha = coralpha * (corstr != "independence"), corstr = corstr,
       ## proximal effect coefficients
       beta0 = setNames(beta0, c("one", "tmod", "base", "state", "lag1a")),
       ## delayed effect coefficients
       beta1 = setNames(beta1, c("one", "lag1tmod", "base", "lag1state")),
       ## treatment probability model coefficients
       eta = setNames(eta, c("one", "base", "state", "lag1a", "lag1y")),
       ## exogenous or time-invariant main effects
       mu = setNames(mu, c("one", "ty", "base")),
       ## time-varying main effects, centered and proximal
       theta0 = setNames(theta0, c("avail", "state")),
       ## time-varying main effects, centered and delayed
       theta1 = setNames(theta1, c("lag1avail", "lag1state")),
       ## availability model coefficients
       coef.avail = setNames(coef.avail, c("one", "tavail", "lag1a", "lag1y")),
       ## binary state model coefficients
       coef.state = setNames(coef.state,
                            c("one", "tstate", "base", "lag1state", "lag1a")),
       ## functions of time in the main effect, proximal effect,
       ## availability model, and binary state model
       tfun = setNames(tfun, c("ty", "tmod", "tavail", "tstate")))
}
