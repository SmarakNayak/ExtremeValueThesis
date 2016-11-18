require(stabledist)
require(fExtremes)
require(plyr)
# assume Mittag-Leffler waiting times
shape <- 0.99
scale <- 1
n = 50000
# generate stables: 
y <- rstable(n = n, alpha = shape, beta = 1, 
             gamma = (cos(shape * pi / 2))^(1/shape), delta = 0, pm = 1)
# exponentials:
x <- rexp(n = n, rate = 1)
# Mittag-Lefflers:
WW <- y * x^(1/shape) * scale
# Event times: 
TT <- cumsum(WW)
# magnitudes of events:
JJ = rgev(n, xi = 0.3, mu = 0, beta = 1)
kMin <- 5
kMax <- min(10000, n)
idxJ <- order(JJ, decreasing = T)[1:kMax]
source(file = "MittagLefflerEstimation.R", local = TRUE)
source(file = "extractDurations.R", local = TRUE)
estimates <- ldply(.data = seq(kMin,kMax), function(k){
  est <- ml.par.est.delta(Tell(TT,idxJ,k),0.05)
  return(c(est$nu, est$CInu, est$delta, est$CIdelta, k))
})
# beta = shape; sigma = scale; topk = number of values used in estimate
names(estimates) <- c("beta","betaL","betaH","delta","deltaL","deltaH","k")

# plots
par(mfrow=c(2,1))

# plot estimates of tail parameter beta
plot(estimates$k,estimates$beta, type="l",ylab= "beta", xlab = "k", ylim = c(0,1.1), 
     main="ML shape parameter")
lines(estimates$k,estimates$betaH, type="l", lty =2)
lines(estimates$k,estimates$betaL, type="l", lty =2)
abline(h = shape, lty = 3)

# plot transformed estimates of scale parameter delta
v <- n/estimates$k
trans.delta <- estimates$delta / v^(1/shape)
trans.delta1 <- estimates$delta / v^(1/(shape+0.02))
trans.delta2 <- estimates$delta / v^(1/(shape-0.02))

plot(estimates$k, trans.delta, type="l", ylab = "delta", xlab = "k",
     main = "ML scale parameter", ylim = range(c(trans.delta1, trans.delta2)))
lines(estimates$k, trans.delta1, lty = 2)
lines(estimates$k, trans.delta2, lty = 2)
abline(h = scale, lty = 3)
