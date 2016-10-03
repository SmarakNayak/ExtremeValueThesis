## Visualisation

Vipaka <- read.csv("Data/Vipaka.csv")
Vipaka <- Vipaka[28:dim(Vipaka)[1], ]
require(parsedate)
require(POT)
parse_iso_8601(Vipaka$time) -> Vipaka$time
n=dim(Vipaka)[1]-27
plot(Vipaka$time, Vipaka$mag, type='h', ylim=c(0,9))

threshold <- function(a){
  Vipaka[Vipaka$mag >= a, ]
}

par(mfrow=c(3,2))
for(a in c(4.5, 5, 5.5,6)){
  potData <- threshold(a)
  plot(potData$time, potData$mag, 'h', ylim = c(0,9))
}

# ML and GP fit
#source("stabilityPlotsV2.R")
TT=Vipaka$time
JJ=Vipaka$mag

# this cutoff translates to this many magnitudes:
epsMax <- 0.2
m <- ceiling(epsMax * n)
idxJ <- order(JJ, decreasing = T)[1:m]

Tell <- function(TT,idxJ,m){
  m=ceiling(m)
  thinT <- sort(TT[idxJ[1:m]])
  diff(c(as.numeric(parse_iso_8601("1970-01-26T10:01:20.000Z")),thinT))
}

source("MittagLefflerDeltaEstimation.R")
# creates a dataframe with point estimates and confidence intervals
# of the Mittag-Leffler parameters mu and nu:
estimates <- ldply(.data = seq(25,m), function(k){
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

# eps := fraction of magnitudes above threshold
estimates$eps <- estimates$topk / n

#plot estimates of scale parameter delta
plot(estimates$topk,estimates$delta, type="l",ylab= "delta", xlab = "k", 
     ylim = c(0,estimates$deltaH[1]), main="ML scale parameter")
lines(estimates$topk,estimates$deltaH, type="l", lty =2)
lines(estimates$topk,estimates$deltaL, type="l", lty =2)

## plot with estimated beta:
#plot(estimates$eps, estimates$delta * (-log(1-estimates$eps))^(1/estimates$beta), type="l", ylim=c(0,10), xlab = "eps", ylab = "b(n) (beta unknown)", main="b(n)")

plot(estimates$eps, estimates$b, type="l", ylim=c(0,max(estimates$b)), xlab = "eps", ylab = "b(n) (beta unknown)", main="b(n)")
lines(estimates$eps,estimates$bH, type="l", lty =2)
lines(estimates$eps,estimates$bL, type="l", lty =2)

## plot with fixed beta estimate:
beta=estimates$beta[m-25]

estimates$bKnown<-estimates$delta * (-log(1-estimates$eps))^(1/beta)
estimates$bKnownH<-estimates$bKnown+(estimates$deltaH-estimates$delta)*(-log(1-estimates$eps))^(1/beta)
estimates$bKnownL<-estimates$bKnown-(estimates$delta-estimates$deltaL)*(-log(1-estimates$eps))^(1/beta)

plot(estimates$eps, estimates$bKnown, type="l", ylim=c(0,max(estimates$bKnownH)), xlab = "eps", ylab = "b(n) (beta constant)", main="b(n)")
lines(estimates$eps,estimates$bKnownH, type="l", lty =2)
lines(estimates$eps,estimates$bKnownL, type="l", lty =2)

#Generalised Pareto Estimate

GPestimates <- ldply(.data = seq(25,m), function(k){
  l=JJ[idxJ[k]]
  est <- fitgpd(JJ[idxJ],l,est = "mle")
  scaleCI=gpd.fiscale(est,0.95)
  shapeCI=gpd.fishape(est,0.95)
  return(c(est$fitted.values[[1]], scaleCI, est$fitted.values[[2]], shapeCI, k))
})
names(GPestimates) <- c("scaleEst","scaleL","scaleH","shapeEst","shapeL","shapeH","topk")

plot(GPestimates$topk,GPestimates$scaleEst, type="l",ylab= "sigma", xlab = "k", 
     ylim = c(0,2), main="GP scale parameter")
lines(GPestimates$topk,GPestimates$scaleH, type="l", lty =2)
lines(GPestimates$topk,GPestimates$scaleL, type="l", lty =2)

plot(GPestimates$topk,GPestimates$shapeEst, type="l",ylab= "xi", xlab = "k", 
     ylim = c(-1,1), main="GP shape parameter")
lines(GPestimates$topk,GPestimates$shapeH, type="l", lty =2)
lines(GPestimates$topk,GPestimates$shapeL, type="l", lty =2)
