## demonstrate models for proximal and delayed treatment effects

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

## fit with 'geeglm'
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

## fit with 'glm'
## nb: in the 'data' argument, a data frame containing a subject identifer must
##     be provided (although it need not be named 'id')
system.time(fitpd.glm <- glm(a ~ lag1a * state, weights = avail,
                             family = "binomial", data = d, subset = time > 3))

## make 'glm' output more like that of 'geeglm'
## nb: this step is necessary for variance estimation later on
fitpd.glm <- glm2gee(fitpd.glm, id)
## nb: consider only the coefficients, as this fit ignores repeated measures
fitpd.glm$coefficients

## --- treatment probability model for weight numerator

## fit with 'geeglm'
fitpn <- geeglm(a ~ 1, id = id, weights = avail, family = "binomial", data = d,
                subset = time > 2, scale.fix = TRUE)
summary(fitpn)

## fit with 'glm'
fitpn.glm <- glm(a ~ 1, weights = avail, family = "binomial", data = d,
                 subset = time > 2)
fitpn.glm <- glm2gee(fitpn.glm, id)
fitpn.glm$coefficients

## --- calculate weights

## nb: fitted values from 'geeglm' or 'glm' can be used interchangeably here
d$pd <- ifelse(d$avail == 0, 0, ifelse(d$time > 3, fitpd$fitted.values, NA))
d$pn <- ifelse(d$avail == 0, 0, ifelse(d$time > 3, fitpn$fitted.values, NA))
d$w <- with(d, ifelse(avail == 0, 0, ifelse(a == 1, pn/pd, (1 - pn)/(1 - pd))))
d$lag1pd <- with(d, delay(id, time, pd))
d$lag1pn <- with(d, delay(id, time, pn))
d$lag1w <- with(d, delay(id, time, w))

## --- estimate the proximal treatment effect

## fit with 'geeglm'
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
##     GENMOD)
estimate(fit1)

## fit with 'lm'
fit1.lm <- lm(y ~ I(time%%2) + varstate + lag1a + state * I(a - pn),
              weights = w, data = d, subset = time > 4)
fit1.lm <- glm2gee(fit1.lm, id)
fit1.lm$vcov <- vcov(fit1.lm, pn = fitpn.glm, pd = fitpd.glm,
                     label = "I(a - pn)")
estimate(fit1.lm)

## --- estimate the delayed treatment effect

## fit with 'geeglm'
fit2 <- geeglm(y ~ I(time%%2) + lag1state + lag1varstate + I(lag1a - lag1pn),
               weights = lag1w, id = id, data = d, subset = time > 5,
               scale.fix = TRUE)
fit2$vcov <- vcov(fit2, pn = fitpn, pd = fitpd, label = "I(lag1a - lag1pn)")
estimate(fit2)

## fit with 'lm'
fit2.lm <- lm(y ~ I(time%%2) + lag1state + lag1varstate + I(lag1a - lag1pn),
              weights = lag1w, data = d, subset = time > 5)
fit2.lm <- glm2gee(fit2.lm, id)
fit2.lm$vcov <- vcov(fit2.lm, pn = fitpn, pd = fitpd,
                     label = "I(lag1a - lag1pn)")
estimate(fit2.lm)
