## demonstrate models for proximal and delayed treatment effects using 'geeglm',
## the model fitting function from the contributed 'geepack' R package

## load functions needed to generate some data
system("R CMD SHLIB rsnmm.c")
dyn.load(if (Sys.info()["sysname"] == "Windows") "rsnmm.dll" else "rsnmm.so")
library("zoo")
source("xzoo.R")
source("rsnmm.R")

set.seed(0)
d <- rsnmm(n = 50, tmax = 200, beta1 = c(-0.1, 0, 0, 0),
           coef.avail = c(log(9), 0, 0, 0))

## define extra variables, using functions from xzoo.R:
## - variation among current and up to the past 2 states
d$varstate <- with(d, roll(id, time, state, width = 3, FUN = var))
## - variation up to the past 3 states
d$lag1varstate <- with(d, delay(id, time, varstate))

## nb: for a given row in 'd'...
##     'time' indexes the treatment occasion
##     'a' is the corresponding treatment indicator
##     'y' is the corresponding proximal response
##     (this is the same format often used for longitudinal data)
d <- subset(d, time > 0)
head(d)

## load functions needed for variance estimation
library("geepack")
source("xgeepack.R")

## --- treatment model (for the weight denominator)

## nb: - weight by availability status so that the variance estimation functions
##       can easily recover the estimating function
##     - omit earlier observations from the data to avoid using initial values
##       in the lagged treatment status ('lag1a')
##     - 'geeglm' approach will be very slow for a large number of subjects and
##       treatment occasions; in this case use 'glm'/'lm' as follows...
system.time(fitpd <- geeglm(a ~ lag1a * state, id = id, weights = avail,
                            family = "binomial", data = d, subset = time > 3,
                            scale.fix = TRUE))
summary(fitpd)

## nb: in the 'data' argument, a data frame containing a subject identifer must
##     be provided (although it need not be named 'id')
system.time(fitpd.glm <- glm(a ~ lag1a * state, weights = avail,
                             family = "binomial", data = d, subset = time > 3))

## --- treatment probability model for weight numerator

fitpn <- geeglm(a ~ 1, id = id, weights = avail, family = "binomial", data = d,
                subset = time > 2, scale.fix = TRUE)
summary(fitpn)

## --- calculate weights

d$pd <- d$pn <- NA
d[rownames(fitpd$fitted.values), "pd"] <- fitpd$fitted.values
d[rownames(fitpn$fitted.values), "pn"] <- fitpn$fitted.values
d[d$avail == 0, c("pd", "pn")] <- 0
d$w <- with(d, ifelse(avail == 0, 0, ifelse(a == 1, pn/pd, (1 - pn)/(1 - pd))))
d$lag1pd <- with(d, delay(id, time, pd))
d$lag1pn <- with(d, delay(id, time, pn))
d$lag1w <- with(d, delay(id, time, w))

## --- estimate the proximal treatment effect

## nb: - any moderators (like 'state') must be included in the regression model
##       via the '*' or ':' *formula* operators
##     - omit earlier observations from the data to avoid basing the state
##       variance variable on initial values or too few values
##     - any observations used to fit the treatment probability model(s), but not
##       the proximal response model should correspond to *earlier* treatment
##       occasions
fit1 <- geeglm(y ~ I(time%%2) + varstate + lag1a + state * I(a - pn),
               id = id, weights = w, data = d, subset = time > 4,
               scale.fix = TRUE)

## adjust variance estimates for estimation of treatment probabilities
## nb: - depending on the 'pn' and 'pd' arguments specified, 'vcov' can handle
##       any combination of centering and weighting
##     - here the 'label' argument is the term label corresponding to the main
##       treatment effect
fit1$vcov <- vcov(fit1, pn = fitpn, pd = fitpd, label = "I(a - pn)")

## summarize the model fit
## nb: 'estimate' can more generally consider linear combinations of regression
##     coefficients, similar to the CONTRAST or ESTIMATE statements in SAS PROC
##     GENMOD
estimate(fit1)
estimate(fit1, rbind("Proximal effect in state -1" = c(rep(0, 5), 1, -1),
                     "Proximal in state 1" = c(rep(0, 5), 1, 1)))

## --- estimate the delayed treatment effect

fit2 <- geeglm(y ~ I(time%%2) + lag1state + lag1varstate + I(lag1a - lag1pn),
               weights = lag1w, id = id, data = d, subset = time > 5,
               scale.fix = TRUE)
fit2$vcov <- vcov(fit2, pn = fitpn, pd = fitpd, label = "I(lag1a - lag1pn)")
estimate(fit2)
