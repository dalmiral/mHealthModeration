############
## README ## 
############
## This code implements the weighted and centered estimator in Boruvka, Almirall, Witkiewitz and Murphy.
## This code is useful when the data arises from a micro-randomized trial in which the randomization
## probabilities are known. In the paper, the known randomization probabilities are p( A_t | H_t ); they are used in 
## the denominator of the weights. In this code, the numerator of the weights, \tilde{p}(A_t | S_t) are estimated
## using logistic regression, and this code adjusts the standard errors of \beta for the fact that \tilde{p} is
## is estimated. Recall that S_t is a subset of H_t.
############
## The MRT data should be in person-period format (long data format) of size n*tmax. 
## There are n subjects, with tmax time points (note: all n subjects have data at 1:tmax time points).
## "a" is treatment; "y" is outcome; "id" is the unique person identifier.
## For the analysis of any new MRT data set, the analyst might change the following lines:
  ## 24, set your working directory where your data is stored
  ## 34, load new data set from working directory
  ## 44-45, the known randomization probability
  ## 49, design matrix for estimating numerator probability, \tilde p
  ## 56, to incorporate availability (I_t in paper), multiply the weights by the availability indicator
  ## 60-62, design matrix for yt+1 (q-size working model + p-size causal model)

## clean all objects in R and set working directory
rm(list=ls(all=TRUE))
setwd("/Users/dalmiral/Dropbox")

## need an expit calculator
expit <- function(x, eta)
{
  e <- as.vector(exp(x %*% eta))
  e / (1 + e)
}

## load MRT data and get sample size (n) and # of time points (tmax)
d <- read.table(file="fakeMRTdata.txt",sep=",")  
tmax <- length( unique(d$time) )
n <- length( unique(d$id) )
small <- if (n <= 50) {TRUE} else {FALSE} ## do the Mancl and DeRouen small sample adjustment to SEs?

  ##############################################
  ## weighting-and-centering estimator #########
  ##############################################
  ## set the randomization probabilities, the "denominator" probabilities
  ## here, the randomization probs are known
  d$denprob <- .5                                             ## p( 1 | H_t ) -- this is the known randomizaiton probability
  d$den     <- d$a*( d$denprob ) + (1-d$a)*( 1-d$denprob )    ## p( A_t | H_t ) -- this is the known one for At=0 and At=1 people
  ##################################
  ## get numerator and center action, the "numerator" probabilities
  ## these are estimated here with ndesign "design matrix"
  ndesign <- cbind(Int = rep(1,n*(tmax) ) )                   ## in this example, it is a intercept-only model
  glmn <- glm.fit(x = ndesign, y = d$a, weights = rep(1,n*(tmax) ), family = binomial() )
  d$numprob <- expit( x = ndesign, eta = glmn$coefficients )  ## \tilde{p}(1 | S_t)
  d$num     <- d$a*( d$numprob ) + (1-d$a)*( 1-d$numprob )    ## \tilde{p}(A_t | S_t)
  d$ac      <- d$a - d$numprob                                ## A_t - \tilde{p}(1 | S_t)
  ##################################
  ## define the weights
  d$wt      <- d$num / d$den    ## multiply this by availability indicator, if appropriate for your application
  ##################################
  ## get weighted estimator
  ## design matrix for y is of size n*tmax x (qq+pp)
  ydesign <- cbind(Int = 1, St=d$St, Ac=d$ac) 
  qq <- 2                                         ## this is q (dimention of alpha -- working model)
  pp <- 1                                         ## this is p (dimension of beta -- causal model)
  if( (qq+pp)!=dim(ydesign)[2] ) stop("q+p must match the dimension of ydesign.")
  lm.wt <- lm.wfit(x = ydesign, y = d$y, w = d$wt)
  ##################################################
  ## get standard errors for weighting-and-centering
  #### bread ####
  Uw <- mapply( function(DD,ww,rr)   t(DD) %*% diag(ww) %*% rr,  
                DD       = split.data.frame( ydesign, d$id ),
                ww       = split(as.vector( d$wt ), d$id),
                rr       = split(d$y - lm.wt$fitted.values, d$id),
                SIMPLIFY = FALSE
  )
  Uw <- do.call("rbind", lapply(Uw, FUN=t) )
  dUw <- mapply(function(DD, ww) t(DD) %*% diag(ww) %*% DD,
                DD       = split.data.frame( ydesign, d$id ),
                ww       = split(as.vector( d$wt ), d$id),
                SIMPLIFY = FALSE
  )
  sumdotUw <- Reduce("+", dUw) ## sum of the derivative of the y estimating function
  EdUw <- sumdotUw / n         ## average of the derivative of the y estimating function
  b <- solve(EdUw)
  #### meat ####
  Un <- mapply( function(DD,rr)   t(DD) %*% rr,    ## numerator estimating function
                DD       = split.data.frame( ndesign, d$id ),
                rr       = split(d$a - d$numprob, d$id),
                SIMPLIFY = FALSE
  )
  Un <- do.call("rbind", lapply(Un, FUN=t) )
  dUn <- mapply(function(DD, ww, prob) t(DD) %*% diag( prob*(1-prob) ) %*% DD, ## derivative of the numerator estimating function
                  DD       = split.data.frame( ndesign, d$id ),
                  prob     = split(d$numprob, d$id),
                  SIMPLIFY = FALSE
  )
  sumdotUn <- Reduce("+", dUn) ## sum of the derivative of the numerator estimating function
  EdUn <- sumdotUn / n         ## average of the derivative of the numerator estimating function
  EdUn.inv <- solve(EdUn)
  Hii <- mapply( function(DDy, ww) DDy %*% solve(sumdotUw) %*% t(DDy) %*% diag(ww), ## tmax x tmax hat matrix for each person
                 DDy      = split.data.frame( ydesign, d$id ),
                 ww       = split(as.vector( d$wt ), d$id),
                 SIMPLIFY = FALSE
  )
  if (small) {  ## small = TRUE, then replace the residual with an adjustment using the hat matrix
    tilderesidualy <- mapply( function(hatmat,rry) solve(diag(1,tmax) - hatmat) %*% rry , ## new residual for each person (scaled by the hat matrix)
                   hatmat   = Hii,
                   rry      = split(d$y - lm.wt$fitted.values, d$id),
                   SIMPLIFY = FALSE
    )
    ## re-define the Uw with the adjustment
    Uw <- mapply( function(DD,ww,rr)   t(DD) %*% diag(ww) %*% rr,  
                  DD       = split.data.frame( ydesign, d$id ),
                  ww       = split(as.vector( d$wt ), d$id),
                  rr       = tilderesidualy,
                  SIMPLIFY = FALSE
    )
    Uw <- do.call("rbind", lapply(Uw, FUN=t) )
  } else {
    tilderesidualy <- split(d$y - lm.wt$fitted.values, d$id)
    ## leave Uw as is.
  }
  Sigwn1 <- mapply( function(DDy,DDn,ww,rry,rrn)   t(DDy) %*% diag(ww) %*% diag( as.vector(rry) ) %*% ( diag(rrn) %*% DDn ),
                   DDy      = split.data.frame( ydesign, d$id ),
                   DDn      = split.data.frame( ndesign, d$id ),
                   ww       = split(as.vector( d$wt ), d$id),
                   rry      = tilderesidualy,
                   rrn      = split(d$a - d$numprob, d$id),
                   SIMPLIFY = FALSE
  )
  Sigwn2 <- mapply( function(DDy,DDn,ww,rry,rrn)   t(DDy) %*% diag(ww) %*% diag( as.vector(rry) ) %*%  diag(rrn) %*% DDn ,
                    DDy      = split.data.frame( cbind( matrix(0,n*(tmax),qq), -1*d$numprob*ydesign[,(qq+1):(qq+pp)]/(d$a - d$numprob) ), d$id ),
                    ww       = split(as.vector( d$wt ), d$id),
                    rry      = tilderesidualy,
                    rrn      = split(1 - d$numprob, d$id),
                    DDn      = split.data.frame( ndesign, d$id ),
                    SIMPLIFY = FALSE
  )
  if(pp==1) { ## hack to get just the causal part without centered action
    justcausal <- ydesign[,(qq+1):(qq+pp)]  * lm.wt$coefficients[(qq+1):(qq+pp)] / (d$a - d$numprob)
  } else { 
    justcausal <- (ydesign[,(qq+1):(qq+pp)] %*% lm.wt$coefficients[(qq+1):(qq+pp)] ) / (d$a - d$numprob)
  }
  Sigwn3 <- mapply( function(DDy,rrn,DDn,ptilde,causal,ww)   t(DDy) %*% diag( ptilde ) %*% diag( causal ) %*% diag( ww ) %*%  diag(rrn) %*% DDn ,
                    DDy      = split.data.frame( ydesign, d$id ),
                    rrn      = split(1 - d$numprob, d$id),
                    DDn      = split.data.frame( ndesign, d$id ),
                    ptilde   = split(as.vector(d$numprob), d$id),
                    causal   = split( justcausal, d$id),
                    ww       = split(as.vector( d$wt ), d$id),
                    SIMPLIFY = FALSE
  )
  Sigwn <- Reduce("+",Sigwn1) / n + Reduce("+",Sigwn2) / n + Reduce("+",Sigwn3) / n
  psi <- Uw + Un %*% EdUn.inv %*% t(Sigwn)
  m <- ( t(psi) %*% psi ) / n
  #### sandwich ####
  varcov <- ( b %*% m %*% t(b) ) / n          ## this one adjusts for numerator
  m <- ( t(Uw) %*% Uw ) / n
  varcov.noNDadj <- ( b %*% m %*% t(b) ) / n  ## this one does NOT adjust for numerator
  #############################################################
  ## get 95% confidence intervals and p-value
  ## use the adjusted SE for confidence intervals and p-values
  se.est <- sqrt( diag(varcov) )
  lcl <- lm.wt$coefficients - se.est * qt(.975, df = n - qq - 1)
  ucl <- lm.wt$coefficients + se.est * qt(.975, df = n - qq - 1)
  pvalue <- pf( (lm.wt$coefficients/se.est)^2, lower.tail = FALSE, df1 = 1, df2 = n - qq - 1) 
  #############################################################
  ## print results
  res <- cbind(
    est=lm.wt$coefficients, 
    unadjusted.se=sqrt( diag(varcov.noNDadj) ), ## unadjusted SE reported but not used for confidence intervals or pvalues
    adjusted.se=se.est, 
    lower95CI=lcl, 
    upper95CI=ucl, 
    pval=pvalue
    ) 
  print( round( res, 4) )
  
  


## eof
