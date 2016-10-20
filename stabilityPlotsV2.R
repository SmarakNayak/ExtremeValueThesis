require(stabledist)
require(fExtremes)
require(plyr)
require(POT) ##install.packages("POT", repos="http://R-Forge.R-project.org") 
# number of observations
n = 100000
# tail parameter of waiting times
beta = 0.4
# scaling constant for "unit" stable under parametrisation pm = 1 below
sigma = (cos(beta * pi / 2))^(1 / beta)
# norming sequence
b.n=n^(-1/beta)
# times of events:
TT = cumsum(rstable(n = n, alpha = beta, beta = 1, 
                    gamma = sigma, delta = 0, pm = 1)) * b.n
# magnitudes of events (distribution irrelevant)
JJ = rgev(n, xi = 0.3, mu = 0, beta = 1)
#Restrict attention to unit interval
# JJ <- JJ[TT < 1]
# TT <- TT[TT < 1]

## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs
source("MittagLefflerEstimation.R")
source("MittagLefflerDeltaEstimation.R")

# consider the cutoff at the top epsMax values:
epsMax <- 0.01
# this cutoff translates to this many magnitudes:
m <- ceiling(epsMax * n)
# indices of the largest jumps:
idxJ <- order(JJ, decreasing = T)[1:m]
source("extractDurations.R")
# creates a dataframe with point estimates and confidence intervals
# of the Mittag-Leffler parameters mu and nu:
estimates <- ldply(.data = seq(50,m), function(k){
  est <- ml.par.est.delta(Tell(TT,idxJ,k),0.05,1-k/n)
  return(c(est$nu, est$CInu, est$delta, est$CIdelta, est$b,est$CIb, k))
})
# beta = shape; sigma = scale; topk = number of values used in estimate
names(estimates) <- c("beta","betaL","betaH","delta","deltaL","deltaH","b","bL","bH" ,"topk")

par(mfrow=c(2,3))
# plot estimates of tail parameter beta
plot(estimates$topk,estimates$beta, type="l",ylab= "beta", xlab = "k", ylim = c(0,1), main="ML tail parameter")
lines(estimates$topk,estimates$betaH, type="l", lty =2)
lines(estimates$topk,estimates$betaL, type="l", lty =2)
abline(h = beta, lty = 3)

# eps := fraction of magnitudes above threshold
estimates$eps <- estimates$topk / n
estimates$truedelta <- (-log(1-estimates$eps))^-(1/beta)*b.n
#plot estimates of scale parameter delta
plot(estimates$topk,estimates$delta, type="l",ylab= "delta", xlab = "k", 
     main="ML scale parameter")
lines(estimates$topk,estimates$deltaH, type="l", lty =2)
lines(estimates$topk,estimates$deltaL, type="l", lty =2)
lines(estimates$topk,estimates$truedelta,type="l",lty=3)



## plot with known beta:
estimates$bKnown<-estimates$delta * (-log(1-estimates$eps))^(1/beta)
estimates$bKnownH<-estimates$bKnown+(estimates$deltaH-estimates$delta)*(-log(1-estimates$eps))^(1/beta)
estimates$bKnownL<-estimates$bKnown-(estimates$delta-estimates$deltaL)*(-log(1-estimates$eps))^(1/beta)

plot(estimates$eps, estimates$bKnown, type="l", ylim=c(0,2*b.n), xlab = "eps", ylab = "b(n) (beta known)", main="b(n)")
lines(estimates$eps,estimates$bKnownH, type="l", lty =2)
lines(estimates$eps,estimates$bKnownL, type="l", lty =2)
abline(h = b.n, lty = 3)

## plot with estimated beta:
#plot(estimates$eps, estimates$delta * (-log(1-estimates$eps))^(1/estimates$beta), type="l", ylim=c(0,10), xlab = "eps", ylab = "b(n) (beta unknown)", main="b(n)")

plot(estimates$eps, estimates$b, type="l", ylim=c(0,2*b.n), xlab = "eps", ylab = "b(n) (beta unknown)", main="b(n)")
lines(estimates$eps,estimates$bH, type="l", lty =2)
lines(estimates$eps,estimates$bL, type="l", lty =2)
abline(h = b.n, lty = 3)

#Generalised Pareto Estimate
#GPmleEst <- fitgpd(JJ[idxJ], l, est = "mle")
#mom <- fitgpd(x, 1, est = "moments")
#pwmb <- fitgpd(x, 1, est = "pwmb")
#pwmu <- fitgpd(x, 1, est = "pwmu")
#gpd.fiscale(GPmleEst, conf = 0.95)
#gpd.fishape(GPmleEst,conf=0.95)
#gpd.pfscale(GPmleEst, conf = 0.95)
#gpd.pfshape(GPmleEst,conf=0.95)

GPestimates <- ldply(.data = seq(50,m), function(k){
  l=JJ[idxJ[k]]
  theoreticalScale=1+0.3*(l)
  est <- fitgpd(JJ[idxJ],l,est = "mle")
  scaleCI=gpd.fiscale(est,0.95)
  shapeCI=gpd.fishape(est,0.95)
  return(c(est$fitted.values[[1]], scaleCI,theoreticalScale, est$fitted.values[[2]], shapeCI, k))
})
names(GPestimates) <- c("scaleEst","scaleL","scaleH","theoreticalScale","shapeEst","shapeL","shapeH","topk")

plot(GPestimates$topk,GPestimates$scaleEst, type="l",ylab= "sigma", xlab = "k", 
     ylim = c(0,10), main="GP scale parameter")
lines(GPestimates$topk,GPestimates$scaleH, type="l", lty =2)
lines(GPestimates$topk,GPestimates$scaleL, type="l", lty =2)
lines(GPestimates$topk,GPestimates$theoreticalScale, type="l", lty =3)

plot(GPestimates$topk,GPestimates$shapeEst, type="l",ylab= "xi", xlab = "k", 
     ylim = c(0,1), main="GP shape parameter")
lines(GPestimates$topk,GPestimates$shapeH, type="l", lty =2)
lines(GPestimates$topk,GPestimates$shapeL, type="l", lty =2)
abline(h = 0.3, lty = 3)
