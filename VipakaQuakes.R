## Visualisation

Vipaka <- read.csv("Data/Vipaka.csv")
Vipaka <- Vipaka[-(1:27), ]
require(parsedate)
require(POT)
parse_iso_8601(Vipaka$time) -> Vipaka$time
n=dim(Vipaka)[1]
plot(Vipaka$time, Vipaka$mag, type='h', ylim=c(0,9))

threshold <- function(a){
  Vipaka[Vipaka$mag >= a, ]
}

par(mfrow=c(2,2))
# for(a in c(4.5, 5, 5.5,6)){
#   potData <- threshold(a)
#   plot(potData$time, potData$mag, 'h', ylim = c(0,9))
# }

# ML and GP fit
#source("stabilityPlotsV2.R")
TT=as.vector(Vipaka$time)
JJ=Vipaka$mag

# set lowest threshold at this many values
k_max = n
idxJ <- order(JJ, decreasing = T)[1:k_max]
# set highest threshold at this many values
k_min = 5

source("get_durations.R")
source("MittagLefflerDeltaEstimation.R")

# creates a dataframe with point estimates and confidence intervals
# of the Mittag-Leffler parameters mu and nu:
estimates <- ldply(.data = seq(k_min, k_max), function(k){
  est <- ml.par.est.delta(get_durations(TT,idxJ,k),0.05,1-k/n)
  return(c(est$nu, est$CInu, est$delta, est$CIdelta, k))
})
# beta = shape; sigma = scale; k = number of values used in estimate
names(estimates) <- c("beta","betaL","betaH","delta","deltaL","deltaH" ,"k")


# plot estimates of tail parameter beta
plot(estimates$k,estimates$beta, type="l", ylab= "beta", 
     xlab = "k", ylim = c(0,1), main="ML tail")
JJ[idxJ[k_max]]
axis(side = 3, at = c(k_min, k_max), labels = c(JJ[idxJ[k_min]], JJ[idxJ[k_max]]))
lines(estimates$k,estimates$betaH, type="l", lty =2)
lines(estimates$k,estimates$betaL, type="l", lty =2)

# extract beta (not implemented)
hat_beta = 0.4

# p := fraction of magnitudes above threshold
estimates$p <- estimates$k

delta0 <- estimates$delta * estimates$p^hat_beta
delta01 <- estimates$delta * estimates$p^(hat_beta-0.05)
delta02 <- estimates$delta * estimates$p^(hat_beta+0.05)

#plot estimates of scale parameter delta

plot(estimates$k, delta0, type='l', main = "ML scale")
lines(estimates$k, delta01, type='l', lty=2)
lines(estimates$k, delta02, type='l', lty=2)

#Generalised Pareto Estimate

GPestimates <- ldply(.data = seq(k_min, k_max), function(k){
  l=JJ[idxJ[k]]
  est <- fitgpd(JJ[idxJ],l,est = "mle")
  scaleCI=gpd.fiscale(est,0.95)
  shapeCI=gpd.fishape(est,0.95)
  return(c(est$fitted.values[[1]], scaleCI, est$fitted.values[[2]], shapeCI, k))
})
names(GPestimates) <- c("scaleEst","scaleL","scaleH","shapeEst","shapeL","shapeH","k")

plot(GPestimates$k,GPestimates$shapeEst, type="l",ylab= "xi", xlab = "k", 
     ylim = c(-1,1), main="GP shape")
lines(GPestimates$k,GPestimates$shapeH, type="l", lty =2)
lines(GPestimates$k,GPestimates$shapeL, type="l", lty =2)

plot(GPestimates$k,GPestimates$scaleEst, type="l",ylab= "sigma", xlab = "k", 
     ylim = c(0,2), main="GP scale")
lines(GPestimates$k,GPestimates$scaleH, type="l", lty =2)
lines(GPestimates$k,GPestimates$scaleL, type="l", lty =2)

