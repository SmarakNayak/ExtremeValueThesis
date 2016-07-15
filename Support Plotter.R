require(stabledist)
require(fExtremes)
n = 1000
par(mfrow=c(1,1))

for (b in seq(0.1,0.9,0.1)) {
  beta = b
  sumWaitingTimes = cumsum(rstable(n,beta, 1, pm=1))
  
  jumps = rgev(n, xi = 0.3, mu = 0, beta = 1)
  
  # Now threshold the data
  # Take the m largest observations
  m = 50
  indexesOfLargestJumps = order(jumps, decreasing = T)[1:m]
  #plot(sumWaitingTimes[indexesOfLargestJumps],jumps[indexesOfLargestJumps], type='h',xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumps)))
  
  a=jumps[indexesOfLargestJumps[m]]
  
  ##Cheating here to calculate actual F(a) rather than estimating it
  Fa=pgev(a,xi = 0.3, mu = 0, beta = 1)
  f1=function(x,y,beta) {-log(Fa)*dexp(-log(Fa)*exp((x-y)*beta),1)*beta*exp((x-y)*beta)*dstable(exp(y),beta,1,pm=1)*exp(y)}
  f2=function(x,y) {f1(x,y,beta = beta)}
  
  x1=c(-20:80)
  y1=seq(-20,40,2)
  z<-outer(x1,y1,f2)
  
  persp(x1,y1,z,theta = 30, phi = 30, expand = 0.5, col = "lightblue",ticktype="detailed",main=bquote(beta==.(beta)))
}