require(stabledist)
require(fExtremes)
require(plyr)
# number of observations
n = 50000
# tail parameter of waiting times
beta = 0.7
# times of events:
TT = cumsum(rstable(n,beta, 1, gamma=1, delta=0, pm=1))
# magnitudes of events (distribution irrelevant)
JJ = rgev(n, xi = 0.3, mu = 0, beta = 1)


## Estimate tail index of inter-arrival times, where data are thinned
## out at different cutoffs
source("MittagLefflerEstimation.R")

# consider the cutoff at the top epsMax values:
epsMax <- 0.10
# this cutoff translates to this many magnitudes:
m <- ceiling(epsMax * n)
# indices of the largest jumps:
idxJ <- order(JJ, decreasing = T)[1:m]
# extracts (preceding) durations of the m largest exceedances
Tell <- function(TT,idxJ,m){
  m=ceiling(m)
  thinT <- sort(TT[idxJ[1:m]])
  diff(c(0,thinT))
}
# creates a dataframe with point estimates and confidence intervals
# of the Mittag-Leffler parameters mu and nu:
estimates <- ldply(.data = seq(50,m), function(k){
  est <- ml.par.est(Tell(TT,idxJ,k),0.05)
  return(c(est$nu, est$CInu, est$mu, est$CImu, k))
})
# beta = shape; sigma = scale; topk = number of values used in estimate
names(estimates) <- c("beta","betaL","betaH","mu","muL","muH","topk")

par(mfrow=c(1,3))
# plot estimates of tail parameter beta
plot(estimates$topk,estimates$beta, type="l",ylab= "beta", xlab = "k", 
     ylim = c(0,1), main="tail parameter")
lines(estimates$topk,estimates$betaH, type="l", lty =2)
lines(estimates$topk,estimates$betaL, type="l", lty =2)
# plot estimates of scale parameter delta
estimates$delta <- estimates$mu^{-1/estimates$beta}
# eps := fraction of magnitudes above threshold
estimates$eps <- estimates$topk / n
# plot with known beta:
plot(estimates$eps, estimates$delta * (-log(1-estimates$eps))^(1/beta), type="l", ylim=c(0,10), xlab = "eps", ylab = "delta (beta known)", main="scale parameter")
# plot with estimated beta:
plot(estimates$eps, estimates$delta * (-log(1-estimates$eps))^(1/estimates$beta), type="l", ylim=c(0,10), xlab = "eps", ylab = "delta (beta unknown)", main="scale parameter")


