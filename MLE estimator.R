require(stabledist)
n = 1000
beta = 0.4
sumWaitingTimes = cumsum(rstable(n,beta, 1, pm=1))
require(fExtremes)
jumpSizes = rgev(n, xi = 0.3, mu = 0, beta = 1)

# This is what the original data could look like
par(mfrow=c(1,2))
plot(sumWaitingTimes,jumpSizes,type='h', xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumpSizes)))

# Now threshold the data
# Take the m largest observations
m = 50
indexesOfLargestJumps = order(jumpSizes, decreasing = T)[1:m]
plot(sumWaitingTimes[indexesOfLargestJumps],jumpSizes[indexesOfLargestJumps], type='h',xlim = c(0,max(sumWaitingTimes)), ylim=c(0,max(jumpSizes)))

a=jumpSizes[indexesOfLargestJumps[m]]

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
totalLikelihood=vector("numeric",length = 41)
betaVec=seq(0.1,0.9,0.02)

for (j in c(1:41))
{
  for (i in c(1:m))
  {
    betaHat=betaVec[j]
    likelihood[i]=as.numeric(integrate(f=function(y){f1(lnTa[i],y,betaHat)},lower=-10,upper=10)[1])
  }
  totalLikelihood[j]=prod(likelihood)
}
par(mfrow=c(1,1))
plot(betaVec,totalLikelihood,main=bquote("MLE estimate for data generated with "*beta == .(beta)),xlab=expression(hat(beta)))

x1=c(-6:16)
y1=seq(-7,10,0.1)
z<-outer(x1,y1,f2)

par(mfrow=c(1,1))
persp(x1,y1,z,theta = 30, phi = 30, expand = 0.5, col = "lightblue",ticktype="detailed")