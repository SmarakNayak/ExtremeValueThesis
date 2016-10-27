require(stabledist)
require(fExtremes)
require(plyr)
require(POT) ##install.packages("POT", repos="http://R-Forge.R-project.org") 
# number of observations
n = 50000
# tail parameter of waiting times
alpha = 0.4
# scaling constant for "unit" stable under parametrisation pm = 1 below
sigma = (cos(alpha * pi / 2))^(1 / alpha)
# block number; maximum of this many observations results in good approximation
# by GEV distribution
B = 1
# norming sequence
b.n=B^(1/alpha)
# times of events:
TT = cumsum(rstable(n = n, alpha = alpha, beta = 1, gamma = sigma, delta = 0, pm = 1))/b.n
# magnitudes of events (distribution irrelevant)
xi=0.3
sigmaGEV=1
JJ = rgev(n, xi = xi, mu = 0, beta = sigmaGEV)

## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs
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
# alpha = shape; sigma = scale; topk = number of values used in estimate
names(estimates) <- c("alpha","alphaL","alphaH","delta","deltaL","deltaH","b","bL","bH" ,"topk")

par(mfrow=c(1,1))
estimates$eps <- B * estimates$topk / n
# plot estimates of tail parameter alpha
plot(estimates$eps,estimates$alpha, type="l",ylab= "alpha", xlab = "epsilon", ylim = c(0.3,0.5), main="ML tail parameter")
lines(estimates$eps,estimates$alphaH, type="l", lty =2)
lines(estimates$eps,estimates$alphaL, type="l", lty =2)
abline(h = alpha, lty = 3,col="red")

# eps := fraction of magnitudes above threshold
estimates$truedelta <- (-log(1-estimates$eps))^-(1/alpha)*b.n
#plot estimates of scale parameter delta
plot(estimates$eps,estimates$delta, type="l",ylab= "delta", xlab = "epsilon", main="ML scale parameter")
lines(estimates$eps,estimates$deltaH, type="l", lty =2)
lines(estimates$eps,estimates$deltaL, type="l", lty =2)
lines(estimates$eps,estimates$truedelta,type="l",lty=3,col="red")



## plot with known alpha:
estimates$bKnown<-estimates$delta * (-log(1-estimates$eps))^(1/alpha)
estimates$bKnownH<-estimates$bKnown+(estimates$deltaH-estimates$delta)*(-log(1-estimates$eps))^(1/alpha)
estimates$bKnownL<-estimates$bKnown-(estimates$delta-estimates$deltaL)*(-log(1-estimates$eps))^(1/alpha)

plot(estimates$eps, estimates$bKnown, type="l", ylim=c(0,2*b.n), xlab = "epsilon", ylab = "delta* (alpha known)", main="ML scale parameter")
lines(estimates$eps,estimates$bKnownH, type="l", lty =2)
lines(estimates$eps,estimates$bKnownL, type="l", lty =2)
abline(h = b.n, lty = 3,col="red")

## plot with estimated alpha:
#plot(estimates$eps, estimates$delta * (-log(1-estimates$eps))^(1/estimates$alpha), type="l", ylim=c(0,10), xlab = "eps", ylab = "b(n) (alpha unknown)", main="b(n)")

plot(estimates$eps, estimates$b, type="l", ylim=c(0,2*b.n), xlab = "eps", ylab = "b(n) (alpha unknown)", main="b(n)")
lines(estimates$eps,estimates$bH, type="l", lty =2)
lines(estimates$eps,estimates$bL, type="l", lty =2)
abline(h = b.n, lty = 3,col="red")

#Generalised Pareto Estimate

GPestimates <- ldply(.data = seq(50,m), function(k){
  l=JJ[idxJ[k]]
  theoreticalScale=sigmaGEV+xi*(l)
  est <- fitgpd(JJ[idxJ],l,est = "mle")
  sigmaStar=est$fitted.values[[1]]-xi*l
  scaleCI=gpd.fiscale(est,0.95)
  shapeCI=gpd.fishape(est,0.95)
  return(c(est$fitted.values[[1]], scaleCI,theoreticalScale, est$fitted.values[[2]], shapeCI,sigmaStar,k))
})
names(GPestimates) <- c("scaleEst","scaleL","scaleH","theoreticalScale","shapeEst","shapeL","shapeH","sigmaStar","topk")

GPestimates$eps <- B * GPestimates$topk / n

plot(GPestimates$eps,GPestimates$scaleEst, type="l",ylab= "sigma", xlab = "epsilon", 
     ylim = c(0,10), main="GP scale parameter")
lines(GPestimates$eps,GPestimates$scaleH, type="l", lty =2)
lines(GPestimates$eps,GPestimates$scaleL, type="l", lty =2)
lines(GPestimates$eps,GPestimates$theoreticalScale, type="l", lty =3,col="red")

plot(GPestimates$eps,GPestimates$shapeEst, type="l",ylab= "xi", xlab = "epsilon", ylim = c(0,1), main="GP shape parameter")
lines(GPestimates$eps,GPestimates$shapeH, type="l", lty =2)
lines(GPestimates$eps,GPestimates$shapeL, type="l", lty =2)
abline(h = xi, lty = 3, col="red")

plot(GPestimates$eps,GPestimates$sigmaStar,type="l",ylab= "sigma*", xlab = "epsilon", ylim = c(0,3), main="GP scale parameter")
lines(GPestimates$eps,GPestimates$sigmaStar+(GPestimates$scaleH-GPestimates$scaleEst), type="l", lty =2)
lines(GPestimates$eps,GPestimates$sigmaStar-(GPestimates$scaleEst-GPestimates$scaleL), type="l", lty =2)
abline(h = 1, lty = 3, col="red")
