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
  #Calculating actual Ta's
  cumulativeDurations=c(0,sort(sumWaitingTimes[indexesOfLargestJumps]))
  Ta=diff(cumulativeDurations)
  lnTa=log(Ta)
  
  ## Estimate tail index of inter-arrival times, where data are thinned
  ## out at different cutoffs
  f1=function(x,y,beta) {-log(Fa)*dexp(-log(Fa)*exp((x-y)*beta),1)*beta*exp((x-y)*beta)*dstable(exp(y),beta,1,pm=1)*exp(y)}
  f2=function(x,y) {f1(x,y,beta = beta)}
  
  
  likelihood=vector("numeric",length = m)
  
  l=41
  totalLikelihood=vector("numeric",length = l)
  betaHatVec=seq(0.1,0.9,length.out = l)
  
  for (j in c(1:l))
  {
    for (i in c(1:m))
    {
      betaHat=betaHatVec[j]
      likelihood[i]=as.numeric(integrate(f=function(y){f1(lnTa[i],y,betaHat)},lower=-20,upper=40)[1])
    }
    totalLikelihood[j]=prod(likelihood)
  }
  plot(betaHatVec,totalLikelihood,main=bquote("MLE estimate for data generated with "*beta == .(beta)),xlab=expression(hat(beta)),ylab="Likelihood")
}