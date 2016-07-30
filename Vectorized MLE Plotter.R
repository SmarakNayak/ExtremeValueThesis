require(stabledist)
require(fExtremes)
n=1000
betaGen = 0.4
sumWaitingTimes = cumsum(rstable(n,betaGen, 1, pm=1))

jumps = rgev(n, xi = 0.3, mu = 0, betaGen = 1)

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
gxy=function(x,y,beta) {-log(Fa)*dexp(-log(Fa)*exp((x-y)*beta),1)*beta*exp((x-y)*beta)}
g0y=function(y,beta) {gxy(0,y,beta)}

fy=function(y,beta) {dstable(exp(y),beta,1,pm=1)*exp(y)}
dy=0.2
intLB=-20
intUB=40
integrationMeshForY=seq(intLB,intUB,dy)
noPointsInMesh = length(integrationMeshForY)

meshedLnTa=round(lnTa/dy)*dy
extendedIntegrationMeshForXminusY=seq(-intUB+min(meshedLnTa),-intLB+max(meshedLnTa),dy)

l=41
totalLikelihood=vector("numeric",length = l)
betaHatVec=seq(0.1,0.9,length.out = l)

for (j in c(1:l))
{
  beta=betaHatVec[j]
  F=fy(integrationMeshForY,beta)
  allPossibleG=g0y(extendedIntegrationMeshForXminusY,beta)
  
  xmin=min(meshedLnTa)
  xm=max(meshedLnTa)
  # startIndex=(xm-x1)/dy+1
  # endIndex=(xm-x1)/dy+noPointsInMesh 
  # 
  # G1=allPossibleG[startIndex:endIndex]
  G=matrix(,m,noPointsInMesh)
  itr=1
  for (x in meshedLnTa) 
  {
    startIndex=(x-xmin)/dy+1
    endIndex=(x-xmin)/dy+noPointsInMesh 
    G[itr,]=allPossibleG[startIndex:endIndex]
    itr=itr+1
  }
  likelihoodVector=G%*%F
  logLikelihoodVector=log(likelihoodVector)
  totalLikelihood[j]=sum(logLikelihoodVector)
  #or 
  #totLikelihood=apply(G%*%F, 1, +)
}
plot(betaHatVec,totalLikelihood,main=bquote("MLE estimate for data generated with "*beta == .(betaGen)),xlab=expression(hat(beta)),ylab="Likelihood")
