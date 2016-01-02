### geepack extras

### like geese.fit, but dispense with scale estimation, limit correlation
### structures and give output similar to geeglm
geese.glm <-
  function(x, y, id, offset = rep(0, N), soffset = rep(0, N),
           weights = rep(1, N), waves = NULL, zsca = matrix(1, N, 1),
           zcor = NULL, corp = NULL, control = geese.control(...), b = NULL,
           alpha = NULL, gm = NULL, family = gaussian(), mean.link = NULL,
           variance = NULL, cor.link = "identity", sca.link = "identity",
           link.same = TRUE, scale.fix = TRUE, scale.value = 1,
           corstr = c("independence", "ar1", "exchangeable", "fixed"), ...)
{
  corstr <- match.arg(corstr)
  N <- length(id)
  z <- list()
  z$geese <- geese.fit(x = x, y = y, id = id, offset = offset, soffset = soffset,
                       weights = weights, waves = waves, zsca = zsca,
                       zcor = zcor, corp = corp, control, b = b, alpha = alpha,
                       gm = gm, family = family, mean.link = mean.link,
                       variance = variance, cor.link = cor.link,
                       sca.link = sca.link, link.same = link.same,
                       scale.fix = scale.fix, scale.value = scale.value,
                       corstr = corstr, ...)
  z$geese$X <- x
  if (scale.fix) z$geese$gamma <- 1
  z$y <- y
  z$family <- family
  ## second derivative of mean function (mu) wrt linear predictor (eta),
  ## assuming binomial familly with logit link
  z$family$mu.eta2 <- function(eta) (exp(eta) - exp(2 * eta)) / (1 + exp(eta))^3
  z$id <- z$geese$id <- id
  z$offset <- offset
  z$weights <- weights
  z$coefficients <- z$geese$beta
  z$corstr <- corstr
  if (is.null(z$offset))
    z$linear.predictors <- z$geese$X %*% z$geese$beta
  else z$linear.predictors <- z$offset + z$geese$X %*% z$geese$beta
  z$fitted.values <- z$family$linkinv(z$linear.predictors)
  z$residuals <- y - z$fitted.values
  class(z$geese) <- "geese"
  z
}

### like mu.eta, but second derivative
mu.eta2 <- list()
mu.eta2$binomial <- function(eta) (exp(eta) - exp(2 * eta)) / (1 + exp(eta))^3

### evaluate link-based function at estimate
mu.etahat <- function(x, order = 1) {
  fun <- if (order == 1) x$family$mu.eta
         else x$family$mu.eta2
  as.vector(fun(x$linear.predictors))
}

### extract geeglm's estimating function, returned as n by dim(beta) matrix
estfun.geeglm <- function(x, wcovinv = NULL, small = FALSE, ...)
{
  if (is.null(wcovinv)) wcovinv <- working.covariance(x, invert = TRUE)
  n <- length(unique(x$geese$id))
  T <- length(x$y) / n
  ## apply Mancl and DeRouen's (2001) small sample correction
  lev <- if (small) leverage(x, wcovinv)
         else lapply(1:n, function(i) rep(0, T))
  psi <- mapply(function(D, w, V, r, h) t(D) %*% diag(w) %*% V %*% (r / (1 - h)),
                D = split.data.frame(x$geese$X * mu.etahat(x), x$geese$id),
                w = split(as.vector(x$weights), x$geese$id),
                V = wcovinv,
                r = split(x$y - x$fitted.values, x$geese$id),
                h = lev,
                SIMPLIFY = FALSE)
  do.call("rbind", lapply(psi, t))
}

leverage <- function(x, wcovinv = NULL)
{
  if (is.null(wcovinv)) wcovinv <- working.covariance(x, invert = TRUE)
  b <- bread.geeglm(x, wcovinv)
  mapply(function(D, w, V) diag(D %*% b %*% t(D) %*% diag(w) %*% V),
         D = split.data.frame(x$geese$X * mu.etahat(x), x$geese$id),
         w = split(as.vector(x$weights), x$geese$id),
         V = wcovinv,
         SIMPLIFY = FALSE)
}

### extract bread from geeglm's sandwich variance estimator
### nb: positive definite
bread.geeglm <- function(x, wcovinv = NULL, sum = TRUE, invert = TRUE, ...)
{
  if (is.null(wcovinv))
    wcovinv <- working.covariance(x, invert = TRUE)
  b <- mapply(function(D, w, V) t(D) %*% diag(w) %*% V %*% D,
              D = split.data.frame(x$geese$X * mu.etahat(x), x$geese$id),
              w = split(as.vector(x$weights), x$geese$id),
              V = wcovinv,
              SIMPLIFY = FALSE)
  if (sum) b <- Reduce("+", b)
  if (invert) b <- solve(b)
  b
}

### like bdiag, but make results similar to cbind or rbind
block.diag <- function(...)
{
  l <- list(...)
  as.matrix(do.call("bdiag", l[sapply(l, length) != 0]))
}

### extract meat from geeglm's sandwich variance estimator
### 'x' is the lag ('lag' + 1) GEE for treatment effects
### 'g' is the GEE for treatment probability given history
### 'gn' is the GEE for treatment probability given history subset
meat.geeglm <- function(x, g = NULL, gn = NULL, wcovinv = NULL, lag = 0,
                        trtlabel = "", ...)
{
  if (is.null(wcovinv))
    wcovinv <- working.covariance(x, invert = TRUE)
  ## nb: any small sample correction enabled via '...' is applied
  ## only to the original estimating function
  psi <- estfun.geeglm(x, wcovinv = wcovinv, ...)
  if (!is.null(g) | !is.null(gn)) {
    ## presume 'x' is centered estimator if weights either zero or one
    ## and no weight numerator model is provided
    center <- all(unique(x$weights) %in% 0:1) & is.null(gn)
    if (center) {
      if (x$family$link != "identity")
        stop("Centering is limited to the identity link")
      if (trtlabel %in% attributes(x$terms)$term.labels)
        x$trtcoef <-
          c(0, with(attributes(x$terms),
                    factors[which(rownames(factors) == trtlabel), ]))
      if (is.null(x$trtcoef))
        stop("Use 'trtlabel' to indicate the label of the main treatment effect")
      gn <- NULL
    }
    if (!is.null(g)) {
      if (length(unique(x$geese$id)) != length(unique(g$geese$id)))
        stop("Models 'x' and 'g' should be based on the same sample")
      if (!is.function(g$family$mu.eta2))
        stop("Link derivative 'g$family$mu.eta2' is not a valid function")
      if (!is.null(gn)
          & length(unique(g$geese$id)) != length(unique(gn$geese$id)))
        stop("Models 'g' and 'gn' should be based on the same data")
      ## use fitted probabilities for decision points tbeg:tend,
      ## to keep aligned with times in 'x'
      tbeg <- as.vector(table(g$geese$id) - table(x$geese$id)) + 1 - lag
      if (any(tbeg < 1))
        stop("Too few measurements or clusters in model 'g'")
      tend <- as.vector(table(g$geese$id)) - lag
      wcovinv.g <- working.covariance(g, invert = TRUE)
      ## derivative of 'g' estimating function
      b.fit <-
        mapply(function(D, w, V, r, X, b)
          t(D) %*% diag(w) %*% V %*% diag(r) %*% X - b,
          D = split.data.frame(g$geese$X * mu.etahat(g, 2), g$geese$id),
          w = split(as.vector(g$weights), g$geese$id),
          V = wcovinv.g,
          r = split(g$y - g$fitted.values, g$geese$id),
          X = split.data.frame(g$geese$X, g$geese$id),
          b = bread.geeglm(g, wcovinv = wcovinv.g, sum = FALSE, invert = FALSE),
          SIMPLIFY = FALSE)
      psi.fit <- estfun.geeglm(g, wcovinv.g)
      ## factor in derivative of 'x' estimating function wrt 'g' parameters
      d.fit <-
        if (center)
          mapply(function(D, a, p) D / (a - p),
                 D = split.data.frame(g$geese$X * mu.etahat(g), g$geese$id),
                 a = split(g$y, g$geese$id),
                 p = split(g$fitted.values, g$geese$id),
                 SIMPLIFY = FALSE)
        else
          mapply(function(D, a, p) (-1)^a * D / (p^a * (1 - p)^(1 - a)),
                 D = split.data.frame(g$geese$X * mu.etahat(g), g$geese$id),
                 a = split(g$y, g$geese$id),
                 p = split(g$fitted.values, g$geese$id),
                 SIMPLIFY = FALSE)
    }
    else {
      if (length(unique(x$geese$id)) != length(unique(gn$geese$id)))
        stop("Models 'x' and 'gn' should be based on the same sample")
      psi.fit <- numeric(0)
      b.fit <- d.fit <- split(numeric(0), gn$geese$id)
      tbeg <- as.vector(table(gn$geese$id) - table(x$geese$id)) + 1 - lag
      if (any(tbeg < 1))
        stop("Too few measurements or clusters in model 'gn'")
      tend <- as.vector(table(gn$geese$id)) - lag
    }
    if (!is.null(gn)) {
      if (!is.function(gn$family$mu.eta2))
        stop("Link derivative 'gn$family$mu.eta2' is not a valid function")
      wcovinv.gn <- working.covariance(gn, invert = TRUE)
      ## derivative of 'gn' estimating function
      b.fit <-
        mapply(function(D, w, V, r, X, b, bp)
          block.diag(t(D) %*% diag(w) %*% V %*% diag(r) %*% X - b, bp),
          D = split.data.frame(gn$geese$X * mu.etahat(gn, 2), gn$geese$id),
          w = split(as.vector(gn$weights), gn$geese$id),
          V = wcovinv.gn,
          r = split(gn$y - gn$fitted.values, gn$geese$id),
          X = split.data.frame(gn$geese$X, gn$geese$id),
          b = bread.geeglm(gn, wcovinv.gn, sum = FALSE, invert = FALSE),
          bp = b.fit,
          SIMPLIFY = FALSE)
      psi.fit <- cbind(estfun.geeglm(gn, wcovinv.gn), psi.fit)
      ## factor in derivative of 'x' estimating function wrt 'gn' parameters
      d.fit <-
        mapply(function(D, a, q, dq, dp)
          cbind((-1)^(1 - a) * D / (q^a * (1 - q)^(1 - a)), dp),
          D = split.data.frame(gn$geese$X * mu.etahat(gn), gn$geese$id),
          a = split(gn$y, gn$geese$id),
          q = split(gn$fitted.values, gn$geese$id),
          dp = d.fit,
          SIMPLIFY = FALSE)
    }
    b.fit <- solve(Reduce("+", b.fit))
    tind <- mapply(function(i, j) i:j, tbeg, tend, SIMPLIFY = FALSE)
    ## derivative of 'x' estimating function wrt 'g'/'gn' parameters
    dpsi <-
      if (center)
        mapply(function(D, e, w, V, r, j, D.fit)
          t(D - r * E) %*% V %*% diag(w * e) %*% D.fit[j, ],
          D = split.data.frame(x$geese$X * mu.etahat(x), x$geese$id),
          E = split.data.frame(matrix(x$trtcoef, nrow(x$geese$X),
                                      ncol(x$geese$X), byrow = TRUE),
                               x$geese$id),
          e = split(as.vector(x$geese$X %*% (x$trtcoef * x$coefficients)),
                    x$geese$id),
          w = split(as.vector(x$weights), x$geese$id),
          V = wcovinv,
          r = split(x$y - x$fitted.values, x$geese$id),
          j = tind,
          D.fit = d.fit,
          SIMPLIFY = FALSE)
      else
        mapply(function(D, w, V, r, j, D.fit)
          t(D) %*% V %*% diag(w * r) %*% D.fit[j, ],
          D = split.data.frame(x$geese$X * mu.etahat(x), x$geese$id),
          w = split(as.vector(x$weights), x$geese$id),
          V = wcovinv,
          r = split(x$y - x$fitted.values, x$geese$id),
          j = tind,
          D.fit = d.fit,
          SIMPLIFY = FALSE)
    dpsi <- Reduce("+", dpsi)
    ## meat here is the original estimating function, minus the product of
    ## its derivative wrt weight parameters, bread of 'x', and estimating
    ## function of 'x'
    psi <- t(mapply(function(e, e.fit) e - dpsi %*% b.fit %*% e.fit,
                    e = split(psi, 1:nrow(psi)),
                    e.fit = split(psi.fit, 1:nrow(psi.fit))))
  }
  Reduce("+", lapply(split(psi, 1:nrow(psi)), function(row) row %o% row))
}

### extract geeglm's working covariance matrix
working.covariance <- function(x, invert = FALSE, wcor = NULL)
{
  if (is.null(wcor)) wcor <- working.correlation(x)
  wcov <- mapply(function(a, R, phi) phi * diag(sqrt(a)) %*% R %*% diag(sqrt(a)),
                 a = split(x$family$variance(x$fitted.values), x$geese$id),
                 R = sapply(x$geese$clusz, function(k) wcor[1:k, 1:k],
                            simplify = FALSE),
                 phi = 1 / x$geese$gamma,
                 SIMPLIFY = FALSE)
  if (invert)
    wcov <- lapply(wcov, solve)
  wcov
}

### extract geeglm's working correlation matrix
working.correlation <- function(x, ...)
{
  R <- x$working.correlation
  if (is.null(R)) {
    R <- diag(max(x$geese$clusz))
    alpha <- x$geese$alpha
    if (length(alpha))
      R[lower.tri(R) | upper.tri(R)] <- alpha
    if (x$corstr == "ar1")
      R <- R^abs(col(R) - row(R))
  }
  R
}

var.weighted <- function(x, g, gn)
{
  ## FIXME: account for other mean-link functions
  if (g$family$family == "binomial")
    g$family$mu.eta2 <- mu.eta2$binomial
  if (gn$family$family == "binomial")
    gn$family$mu.eta2 <- mu.eta2$binomial
  b <- bread.geeglm(x)
  m <- meat.geeglm(x, g, gn)
  x$var <- b %*% m %*% t(b)
  x$df <- with(x, length(unique(id)) - length(coef))
  x
}

### summarize linear combinations of regression coefficients
estimate <- function(fit, combos = NULL, var = fit$var, ztest = TRUE,
                     df = length(fit$id) - length(fit$coef))
{
  if (is.null(combos)) {
    combos <- diag(length(fit$coef))
    rownames(combos) <- names(fit$coef)
  }
  if (ztest) {
    qfun <- qnorm
    pfun <- function(q) pchisq(q, df = 1, lower.tail = FALSE)
  }
  else {
    qfun <- function(p) qt(p, df = df)
    pfun <- function(q) pf(q, lower.tail = FALSE, df1 = 1, df2 = df)
  }
  est <- combos %*% fit$coef
  if (is.null(var)) var <- fit$geese$vbeta
  se.est <- sqrt(diag(combos %*% var %*% t(combos)))
  lcl <- est - se.est * qfun(0.975)
  ucl <- est + se.est * qfun(0.975)
  pvalue <- pfun((est/se.est)^2)
  out <- cbind(est, lcl, ucl, se.est, pvalue)
  rownames(out) <- rownames(combos)
  colnames(out) <- c("Estimate", "95% LCL", "95% UCL", "SE", "p-value")
  class(out) <- "estimate"
  out
}

print.estimate <- function(object, digits = min(getOption("digits"), 3),
                           signif.stars = TRUE, ...)
{
  printCoefmat(object, digits = digits, dig.tst = digits,
               signif.stars = signif.stars, has.Pvalue = TRUE,
               eps.Pvalue = 1e-4, ...)
}
