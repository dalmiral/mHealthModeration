library("geepack")
library("Matrix")
library("zoo")

source("xgeepack.R")
source("xzoo.R")
system("R CMD SHLIB rsnmm.c")
dyn.load("rsnmm.so")
source("rsnmm.R")

standardize <- function(x, na.rm = FALSE, scale = TRUE)
  (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)^scale

inside.ci <- function(truth, est, se, cp, qfun = qnorm, ...)
{
  half <- se * qfun(1 - (1 - cp)/2, ...)
  truth >= est - half & truth <= est + half  
}

expit <- function(x, eta)
{
  e <- as.vector(exp(x %*% eta))
  e / (1 + e)
}

dexpit <- function(x, eta)
{
  e <- as.vector(exp(x %*% eta))
  d <- e / (1 + e)^2 * x[, c(1, 1 + which(eta[-1] != 0))]
  if (is.null(ncol(d)))
    d <- matrix(d, ncol = 1)
  d
}

boxcox <- function(x, lambda)
{
  if (lambda[1] == 0)
    log(x + lambda[2])
  else
    ((x + lambda[2])^lambda[1] - 1) / lambda[1]
}

boxcox.inv <- function(x, lambda)
{
  if (lambda[1] == 0)
    exp(x) - lambda[2]
  else
    (lambda[1] * x + 1)^(1/lambda[1]) - lambda[2]
}

params <- function(x, keep = "one", omit = "")
  x[which(x != 0 & !(names(x) %in% omit) | names(x) %in% keep)]

### return names as a quoted list
qnames <- function(x, prefix = NULL, suffix = NULL, keep = "one", omit = "")
{
  if (!is.character(x)) x <- names(x)
  x <- x[x %in% keep | !(x %in% omit)]
  k <- length(x)
  if (!is.null(prefix))
    prefix <- rep(paste(prefix, "*", sep = ""), each = k)
  if (!is.null(suffix))
    suffix <- rep(paste(suffix, "*", sep = ""), each = k)
  x <- paste(prefix, x, suffix, sep = "")
  x <- gsub("^one\\*", "", x)
  x <- gsub("\\*one$", "", x)
  sapply(x, function(n) parse(text = n)[[1]], simplify = FALSE)
}

### return list of zeros of the same length
zeros <- function(x) lapply(x, function(i) 0)

### row-wise outer product
outer.row <- function(u, v = NULL)
{
  m <- ncol(u)
  u <- if (is.null(v)) cbind(u, u)
       else cbind(u, v)
  t(apply(u, 1, function(r) as.vector(r[1:m] %o% r[-(1:m)])))
}
